% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_power.R
\name{plot_power}
\alias{plot_power}
\alias{plotpower}
\title{Power curves for different sampling efforts}
\usage{
plot_power(data, n = NULL, m, method = "both")
}
\arguments{
\item{data}{Object of class "ecocbo_beta" that results from \code{\link[=sim_beta]{sim_beta()}}.}

\item{n}{Defaults to NULL, and then the function computes the number of
samples (n) that results in a sampling effort close to 95\% in power. If
provided, said number of samples will be used.}

\item{m}{Site label to be used as basis for the plot.}

\item{method}{The desired plot. Options are "power", "density" or "both".
"power" plots the power curve, "density" plots the density distribution of
pseudoF, and "both" draws both plots one next to the other.}
}
\value{
If the method is "power", then a power curve in which the selected,
or computed, "n" is marked in red. If the method is "density", then a
density plot for the observed pseudoF values and a line marking the value of
pseudoF that marks the significance level indicated in \code{\link[=sim_beta]{sim_beta()}}.
If the method is "both", then a composite with a power curve and a
density plot side by side.

The value of the selected "m", "n" and the corresponding component of variation
are presented in all methods.
}
\description{
\code{plot_power()} can be used to visualize the power of a study as a
function of the sampling effort. The power curve plot shows that the
power of the study increases as the sample size increases, and the density
plot shows the overlapping areas where \eqn{\alpha} and \eqn{\beta} are
significant.
}
\examples{
plot_power(data = epiBetaR, n = 4, m = 2, method = "both")
plot_power(data = epiBetaR, n = NULL, m = 3, method = "power")
plot_power(data = epiBetaR, n = NULL, m = 3, method = "density")
}
\references{
Underwood, A. J. (1997). Experiments in ecology: their logical
design and interpretation using analysis of variance. Cambridge university
press.

Underwood, A. J., & Chapman, M. G. (2003). Power, precaution,
Type II error and sampling design in assessment of environmental impacts.
Journal of Experimental Marine Biology and Ecology, 296(1), 49-70.
}
\seealso{
\code{\link[=sim_beta]{sim_beta()}}
\code{\link[=scompvar]{scompvar()}}
\code{\link[=sim_cbo]{sim_cbo()}}
}
\author{
Edlin Guerra-Castro (\email{edlinguerra@gmail.com}), Arturo Sanchez-Porras
}
