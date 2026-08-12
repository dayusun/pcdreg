# Turn a PanelCount response and a model matrix into the arrays the C++ core
# expects.
#
# The likelihood involves the covariate trajectory evaluated at every pooled
# examination time up to a subject's own last examination, so the user's rows
# (which are intervals of constant covariates) are expanded onto that pooled
# grid.  Each expanded row is one (subject, grid time) pair, that is, one (i, k)
# with Delta_ik = 1.
prepare_panel <- function(y, X) {
  ord <- order(y[, "id"], y[, "tstop"])
  id <- y[ord, "id"]
  tstart <- y[ord, "tstart"]
  tstop <- y[ord, "tstop"]
  exam <- y[ord, "exam"] == 1
  X <- X[ord, , drop = FALSE]

  used <- sort(unique(id))
  subject <- match(id, used)
  n <- length(used)
  labels <- attr(y, "labels")
  labels <- if (is.null(labels)) as.character(used) else labels[used]

  # Tolerance for recognising that one interval ends where the next begins.
  tol <- 1e-8 * max(1, max(abs(tstop)))

  first <- !duplicated(subject)
  if (any(abs(tstart[first]) > tol)) {
    bad <- which(first)[abs(tstart[first]) > tol][1L]
    stop("Follow-up must start at time 0, but subject ", labels[subject[bad]],
         " starts at ", tstart[bad],
         ". Shift the time scale so that 0 is the start of follow-up.",
         call. = FALSE)
  }
  gap <- !first & abs(tstart - c(0, tstop[-length(tstop)])) > tol
  if (any(gap)) {
    bad <- which(gap)[1L]
    stop("Intervals within a subject must be contiguous, but subject ",
         labels[subject[bad]], " has a gap or overlap at time ", tstart[bad],
         ".", call. = FALSE)
  }
  if (anyDuplicated(cbind(subject, tstop))) {
    stop("Each subject may have at most one row ending at any given time.",
         call. = FALSE)
  }
  nexam_by_subject <- tabulate(subject[exam], n)
  if (any(nexam_by_subject == 0L)) {
    stop("Every subject needs at least one examination, but subject ",
         labels[which(nexam_by_subject == 0L)[1L]], " has none.", call. = FALSE)
  }

  # Follow-up ends at the last examination; later intervals say nothing about
  # the model, because the likelihood only runs over j = 1, ..., J_i.
  subject_f <- factor(subject, levels = seq_len(n))
  last_exam <- as.numeric(tapply(tstop[exam], subject_f[exam], max))
  keep <- tstart < last_exam[subject] - tol
  if (!all(keep)) {
    message("Dropping ", sum(!keep), " row(s) beyond the last examination of ",
            "their subject.")
    subject <- subject[keep]
    subject_f <- subject_f[keep]
    tstop <- tstop[keep]
    X <- X[keep, , drop = FALSE]
  }

  times <- sort(unique(y[y[, "exam"] == 1, "tstop"]))
  K <- length(times)

  # One expanded row per grid time still under observation.
  nactive <- findInterval(last_exam, times)
  subj_exp <- rep.int(seq_len(n), nactive)
  grid_exp <- sequence(nactive)

  # Examination rows, ordered to match the panel intervals they close.
  exam_rows <- y[y[, "exam"] == 1, , drop = FALSE]
  exam_subject <- match(exam_rows[, "id"], used)
  exam_ord <- order(exam_subject, exam_rows[, "tstop"])
  exam_subject <- exam_subject[exam_ord]
  exam_grid <- match(exam_rows[exam_ord, "tstop"], times)

  # Mark the expanded rows that sit at one of their own subject's examinations.
  key <- function(i, k) (as.numeric(i) - 1) * K + k
  is_exam <- key(subj_exp, grid_exp) %in% key(exam_subject, exam_grid)

  # Examination interval index j within subject: one more than the number of
  # that subject's earlier examinations.
  before <- cumsum(is_exam) - is_exam
  j_within <- before - rep.int(before[!duplicated(subj_exp)], nactive) + 1L

  J <- tabulate(exam_subject, n)
  panel <- cumsum(c(0L, J))[subj_exp] + j_within

  # The covariate row in force at each grid time, located within each subject.
  rows_by_subject <- split(seq_along(subject), subject_f)
  times_by_subject <- split(times[grid_exp], factor(subj_exp, levels = seq_len(n)))
  row_index <- unlist(Map(
    function(rows, query) {
      rows[findInterval(query, tstop[rows], left.open = TRUE) + 1L]
    },
    rows_by_subject, times_by_subject
  ), use.names = FALSE)
  if (anyNA(row_index)) {
    stop("The covariate history does not cover the whole follow-up of every ",
         "subject.", call. = FALSE)
  }

  list(
    X = t(X[row_index, , drop = FALSE]),
    subj = as.integer(subj_exp) - 1L,
    grid = as.integer(grid_exp) - 1L,
    panel = as.integer(panel) - 1L,
    dN = as.numeric(exam_rows[exam_ord, "count"]),
    panelsubj = as.integer(rep.int(seq_len(n), J)) - 1L,
    n = n,
    K = K,
    times = times,
    labels = labels,
    nexam = nrow(exam_rows)
  )
}
