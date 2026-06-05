# =============================================================================
# Feature-group tidyselect helpers + group resolution / selection
# =============================================================================
# Three layers, kept deliberately separate:
#
#   1. SCHEMA      FEATURE_PATTERNS (regex per real column family) and
#                  SUPERGROUPS (coarse categorisation). Defined in config.R.
#                  One entry per real family; each family in exactly one
#                  supergroup. Stable; changing it is a schema change.
#
#   2. INTENT      Reusable named selections live in GROUP_BUNDLES (config.R).
#                  A bundle is the reusable form of plotting/modelling intent —
#                  the proper home for what the old `_some` groups were trying
#                  to be. One-off intent is the per-call pick/drop arguments.
#                  Neither touches the schema.
#
#   3. RESOLUTION  resolve_selection() normalises (groups, pick, drop) +
#                  bundles into a (group_keys, pick, drop) triple.
#                  expand_groups() is the group-key-only view of that.
#                  select_features() flattens the triple to actual COLUMN
#                  names present in `df`. Plots that need group identity per
#                  column (e.g. the correlation dotplot) consume the triple
#                  directly. Column selection is thus defined in one place.
#
# --- Bundles -----------------------------------------------------------------
# A GROUP_BUNDLES entry is a list with any of:
#   groups : character vector of group / supergroup / OTHER bundle names
#   pick   : named list  group_key -> columns to KEEP from that group
#   drop   : named list  group_key -> columns to REMOVE from that group
# A bare character vector is accepted as shorthand for list(groups = <vec>).
#
# Merge policy when a caller passes pick/drop AND names a bundle carrying
# pick/drop for the SAME group key: the CALLER wins for that key (per-group
# replacement); the bundle's entry applies only where the caller is silent.
# Resolution order is always pick-then-drop, so a caller drop trims whatever
# a bundle pick produced.
#
# Usage:
#   df %>% select(fg("rnafold_zscores"))
#   df %>% select(all_of(select_features(df, groups = "structure")))
#   select_features(df, groups = "nmd_reported")          # a bundle
#   select_features(df, groups = "nmd",
#                   pick = list(nmd = c("nmd_snv_fragile_codon_density_mrna",
#                                       "nmd_alt_stop_codon_density_mrna")))
# =============================================================================


#' Return a tidyselect spec for a named feature group
#'
#' @param group Character, a key of FEATURE_PATTERNS.
#' @return A tidyselect spec usable inside dplyr::select().
#' @export
fg <- function(group) {
  if (!group %in% names(FEATURE_PATTERNS)) {
    stop("Unknown feature group '", group, "'. Known groups: ",
         paste(names(FEATURE_PATTERNS), collapse = ", "))
  }
  tidyselect::matches(FEATURE_PATTERNS[[group]])
}


#' List the columns that match a named feature group in a given dataframe
#'
#' @param df A dataframe.
#' @param group Character, a key of FEATURE_PATTERNS.
#' @return Character vector of column names, in `df` column order.
#' @export
fg_columns <- function(df, group) {
  if (!group %in% names(FEATURE_PATTERNS)) {
    stop("Unknown feature group '", group, "'.")
  }
  grep(FEATURE_PATTERNS[[group]], names(df), value = TRUE)
}


# --- internal: bundle registry accessor --------------------------------------
.group_bundles <- function() {
  if (exists("GROUP_BUNDLES", inherits = TRUE)) GROUP_BUNDLES else list()
}

# --- internal: normalise a bundle entry to list(groups, pick, drop) ----------
.as_bundle <- function(entry) {
  if (is.character(entry)) entry <- list(groups = entry)
  list(
    groups = if (is.null(entry$groups)) character() else entry$groups,
    pick   = if (is.null(entry$pick))   list()       else entry$pick,
    drop   = if (is.null(entry$drop))   list()       else entry$drop
  )
}


#' Resolve a selection (groups + caller pick/drop) into a normalised triple
#'
#' Expands supergroups and bundles, merges any pick/drop the bundles carry
#' with the caller's (caller wins per group key), and returns the pieces every
#' downstream consumer needs. This is the single source of truth for what a
#' selection means; expand_groups() and select_features() are thin views over
#' it.
#'
#' @param groups Character vector of group / supergroup / bundle names, or NULL
#'   for every FEATURE_PATTERNS key.
#' @param pick   Named list: caller's per-group keep-lists.
#' @param drop   Named list: caller's per-group drop-lists.
#' @return list(groups = <FEATURE_PATTERNS keys>, pick = <named list>,
#'   drop = <named list>).
#' @export
resolve_selection <- function(groups = NULL, pick = list(), drop = list()) {
  bundles <- .group_bundles()

  if (is.null(groups)) {
    return(list(groups = names(FEATURE_PATTERNS), pick = pick, drop = drop))
  }

  out_groups <- character()
  bun_pick   <- list()
  bun_drop   <- list()

  # `seen` guards against bundle self / mutual reference.
  walk <- function(tokens, seen) {
    for (g in tokens) {
      if (g %in% names(SUPERGROUPS)) {
        out_groups <<- c(out_groups, SUPERGROUPS[[g]])
      } else if (g %in% names(bundles)) {
        if (g %in% seen) {
          warning("Bundle '", g, "' is self-referential — cycle broken")
          next
        }
        b <- .as_bundle(bundles[[g]])
        # Bundle pick/drop accumulate; later bundles override earlier ones for
        # the same group key (caller still overrides all of them, below).
        bun_pick[names(b$pick)] <<- b$pick
        bun_drop[names(b$drop)] <<- b$drop
        walk(b$groups, c(seen, g))
      } else if (g %in% names(FEATURE_PATTERNS)) {
        out_groups <<- c(out_groups, g)
      } else {
        warning("Unknown group, supergroup or bundle: '", g, "' — skipped")
      }
    }
  }
  walk(groups, character())

  # Caller wins per group key; bundle entries fill the gaps.
  merged_pick <- utils::modifyList(bun_pick, pick)
  merged_drop <- utils::modifyList(bun_drop, drop)

  list(groups = unique(out_groups), pick = merged_pick, drop = merged_drop)
}


#' Expand group / supergroup / bundle names into FEATURE_PATTERNS keys
#'
#' The group-key-only view of resolve_selection(); pick/drop carried by any
#' named bundles are resolved but not returned (use resolve_selection() or
#' select_features() if you need them).
#'
#' @param groups Character vector, or NULL for every FEATURE_PATTERNS key.
#' @return Character vector of FEATURE_PATTERNS keys.
#' @examples
#' expand_groups()                              # every group
#' expand_groups("structure")                   # all structure groups
#' expand_groups(c("structure", "junctions"))   # mixed
#' @export
expand_groups <- function(groups = NULL) {
  resolve_selection(groups)$groups
}


#' Apply pick/drop refinement to one group's columns (shared semantics)
#'
#' pick = allow-list, honouring caller order; drop = remove from the whole
#' group. Order: pick then drop. Used by both select_features() and the
#' correlation dotplot so the two never diverge.
#'
#' @param members Character vector of columns for a group, in df order.
#' @param pick_g  Character vector or NULL.
#' @param drop_g  Character vector or NULL.
#' @return Refined character vector.
#' @export
refine_group_columns <- function(members, pick_g = NULL, drop_g = NULL) {
  if (!is.null(pick_g)) members <- pick_g[pick_g %in% members]   # caller order
  if (!is.null(drop_g)) members <- setdiff(members, drop_g)      # members order
  members
}


#' Resolve a feature selection into an ordered vector of column names
#'
#' The single entry point for "which columns does this analysis use". Accepts
#' groups / supergroups / bundles plus optional per-group pick (allow-list) and
#' drop (remove). See file header for the bundle data model and merge policy.
#'
#' Missing columns are silently skipped (loaders drop missing data — R5).
#' pick/drop names that match no column are reported via message().
#'
#' @param df     Dataframe from build_dataset() / build_all().
#' @param groups Character vector of group / supergroup / bundle names, or NULL.
#' @param pick   Named list: group key -> columns to keep.
#' @param drop   Named list: group key -> columns to drop.
#' @return Ordered, de-duplicated character vector of column names in `df`.
#' @export
select_features <- function(df, groups = NULL,
                            pick = list(), drop = list()) {
  sel      <- resolve_selection(groups, pick, drop)
  expanded <- sel$groups

  unknown_pick <- setdiff(names(sel$pick), expanded)
  unknown_drop <- setdiff(names(sel$drop), expanded)
  if (length(unknown_pick) > 0) {
    message("select_features: pick names not in resolved groups (ignored): ",
            paste(unknown_pick, collapse = ", "))
  }
  if (length(unknown_drop) > 0) {
    message("select_features: drop names not in resolved groups (ignored): ",
            paste(unknown_drop, collapse = ", "))
  }

  cols <- character()
  for (g in expanded) {
    members <- fg_columns(df, g)

    if (!is.null(sel$pick[[g]])) {
      missing <- setdiff(sel$pick[[g]], members)
      if (length(missing) > 0) {
        message("select_features: pick[['", g, "']] columns not found: ",
                paste(missing, collapse = ", "))
      }
    }
    if (!is.null(sel$drop[[g]])) {
      missing <- setdiff(sel$drop[[g]], members)
      if (length(missing) > 0) {
        message("select_features: drop[['", g, "']] columns not found: ",
                paste(missing, collapse = ", "))
      }
    }

    cols <- c(cols, refine_group_columns(members, sel$pick[[g]], sel$drop[[g]]))
  }

  unique(cols)
}
