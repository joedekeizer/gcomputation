\name{plot.gcsurv}
\alias{plot.gcsurv}
\title{

S3 Method for Plotting an 'gcsurv' Object

}
\description{

Provides graphical summaries for an object returned by the function \code{gc_survival}.
This S3 method offers two types of plots: a calibration plot and a survival curve plot.

}
\usage{

\method{plot}{gcsurv}(x, method="calibration", n.groups=5, pro.time=NULL, ...)

}
\arguments{

\item{x}{An object returned by the function \code{gc_survival}.}

\item{method}{A character string specifying the type of plot. Possible values are \code{calibration} (default) or \code{survival}.}

\item{n.groups}{An integer for the number of groups to divide the predicted probabilities into when \code{method="calibration"} is used. The default is 5.}

\item{pro.time}{The prognostic time at which the calibration plot of the survival probabilities is performed. If \code{NULL} (default), the median of the outcome time variable is used.}

\item{...}{Additional graphical parameters that can be passed to the underlying plot function.}

}
\details{

The \code{plot.gcsurv} function is an S3 method for visualizing results from a \code{gcsurv} object, which is generated by the \code{gc_survival} function.

This plotting method supports two distinct types of visualizations:
\itemize{
\item \code{method="survival"}: This plot compares the overall Kaplan-Meier (KM) estimate of the survival curve from the original dataset against the mean of the survival predictions derived from the G-Computation model. The KM curve is typically displayed in blue, while the mean of the G-computation predictions is shown in red. This plot provides a visual assessment of how well the G-computation model's average survival predictions align with the observed survival experience in the raw data.

\item \code{method="calibration"}: This plot assesses the calibration of the Q-model's predictions used within G-Computation created on all of the raw data.
    \itemize{
    \item It takes the predicted survival probabilities at a specific \code{pro.time}. If \code{pro.time} is not specified, it defaults to the median of the outcome time variable.
    \item These predicted probabilities are then divided into \code{n.groups} (defaulting to 5) by their quantiles, creating different risk strata.
    \item For each group, the mean predicted survival probability is plotted on the X axe against the observed KM survival probability for that group at the specified \code{pro.time} on the Y axe.
    \item 95\% confidence intervals for these observed survival probabilities are also displayed as vertical segments.
    \item An identity line (y=x) is included for reference. Perfect calibration is indicated when observed points are along this line, suggesting that the predicted probabilities accurately reflect the observed survival rates across different risk groups.
    }
}


}
\value{No return value for this S3 method.}
\examples{
data(dataPROPHYVAP)

dataPROPHYVAP$DEATH_J60_num <- ifelse(dataPROPHYVAP$DEATH_J60 == "Yes", 1, 0)
dataPROPHYVAP$GROUP_num <- ifelse(dataPROPHYVAP$GROUP == "Placebo", 0, 1)

formula <- formula(Surv(FOLLOWUP_J60, DEATH_J60_num) ~ GROUP_num * (AGE + SEXE +
              BMI + DIABETES))

gc1 <- gc_survival(formula=formula, method="lasso", data=dataPROPHYVAP,
              group="GROUP_num", param.tune=NULL, boot.type="bcv", cv=10,
              boot.number=500, effect="ATE", progress=TRUE , max.time=10,
              boot.tune = TRUE)

plot(gc1)

#Plot the calibration curve
plot(gc1, method="calibration", n.groups=5)

#Plot the survival curves
plot(gc1, method="survival")

}
\keyword{Plotting}
\keyword{G-computation}
\keyword{Survival}