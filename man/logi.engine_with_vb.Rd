\name{logi.engine_with_vb}
\alias{logi.engine_with_vb}
\title{Override spatstat's logi.engine}
\usage{
logi.engine_with_vb(Q, trend = ~1, interaction, ..., covariates = NULL,
  correction = "border", rbord = reach(interaction), covfunargs = list(),
  allcovar = FALSE, vnamebase = c("Interaction", "Interact."),
  vnameprefix = NULL, justQ = FALSE, savecomputed = FALSE,
  precomputed = NULL, VB = FALSE)
}
\description{
Override spatstat's logi.engine
}

