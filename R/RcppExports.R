# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

worldAnalysis_ <- function(wave, frameshift, fs, f0floor, allowed_range) {
    .Call('_WorldR_worldAnalysis_', PACKAGE = 'WorldR', wave, frameshift, fs, f0floor, allowed_range)
}

worldAnalysis_f0 <- function(wave, Rf0, frameshift, fs, f0floor, allowed_range) {
    .Call('_WorldR_worldAnalysis_f0', PACKAGE = 'WorldR', wave, Rf0, frameshift, fs, f0floor, allowed_range)
}

worldF0Estimation_dio <- function(wave, frameshift, fs, f0floor, allowed_range) {
    .Call('_WorldR_worldF0Estimation_dio', PACKAGE = 'WorldR', wave, frameshift, fs, f0floor, allowed_range)
}

worldF0Estimation_harvest <- function(wave, frameshift, fs, f0floor, f0ceil) {
    .Call('_WorldR_worldF0Estimation_harvest', PACKAGE = 'WorldR', wave, frameshift, fs, f0floor, f0ceil)
}

worldSynthesis_ <- function(world) {
    .Call('_WorldR_worldSynthesis_', PACKAGE = 'WorldR', world)
}

