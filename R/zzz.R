.onAttach <- function(libname, pkgname) {
    packageStartupMessage("fraq ", utils::packageVersion("fraq"))
}
