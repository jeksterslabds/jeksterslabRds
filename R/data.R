#' Galton's Height Data
#'
#' @format A data frame with 898 cases and 6 variables:
#' \describe{
#'   \item{family}{The family that the child belongs to, labeled by the numbers from 1 to 204 and 136A.}
#'   \item{father}{The father's height, in inches.}
#'   \item{mother}{The mother's height, in inches.}
#'   \item{gender}{The gender of the child, male (M) or female (F).}
#'   \item{height}{The height of the child, in inches.}
#'   \item{kids}{The number of kids in the family of the child.}
#'   \item{male}{1 if the child is male. 0 if the child is female.}
#'   \item{female}{1 if the child is female. 0 if the child is male.}
#' }
#' @source Francis Galton, 2017, "Galton height data",
#'   \url{https://doi.org/10.7910/DVN/T0HSJ1}, Harvard Dataverse, V1,
#'   UNF:6:2ty+0YgqR2a66FlvjCuPkQ== \[fileUNF\]
"galton"


#' Mothers' and daughters' Heights
#'
#' Heights of mothers in the United Kingdom under the age of 65 and
#'   one of their adult daughters over the age of 18
#'   collected during the period 1893-1898.
#'
#' @format A dataframe with 1375 rows and two variables
#' \describe{
#'   \item{mheight}{Mother's height, in inches.}
#'   \item{dheight}{Daughter's height, in inches.}
#'  }
#' @docType data
#' @usage data(heights)
#' @keywords datasets
#' @references
#'   Pearson, E. S. and Lee, S. (1903). On the laws of inheritance in man. _Biometrika_, 2, 357-–463.
#' @examples
#' data(heights)
"heights"

#' Wages of Workers
#'
#' Wages and related data for 1289 workers.
#'
#' @format A dataframe with 1289 rows and seven variables
#' \describe{
#'   \item{wages}{Hourly wages in dollars.}
#'   \item{gender}{1 for female, 0 for male.}
#'   \item{race}{1 for nonwhite workers and 0 for white workers.}
#'   \item{union}{1 if in union, 0 otherwise.}
#'   \item{edu}{Education in years.}
#'   \item{exp}{Work experience in years.}
#'   \item{age}{Age in years.}
#'  }
#' @docType data
#' @usage data(wages)
#' @keywords datasets
#' @examples
#' data(wages)
"wages"

#' Wages of Youth
#'
#' Wages and educational attainment of youth in the USA
#'
#' @format A dataframe with 540 rows and 16 variables
#' \describe{
#'   \item{age}{Age in 2002.}
#'   \item{school}{Years of schooling (highest grade completed as of 2002).}
#'   \item{female}{1 for female, 0 for male.}
#'   \item{male}{1 for female, 0 for male.}
#'   \item{black}{1 for blacks.}
#'   \item{hisp}{1 for Hispanic; non-black, non-Hispanic is the reference category.}
#'   \item{exp}{Total out-of-school work experience in years as of the 2002 interview.}
#'   \item{married}{1 if married, 0 otherwise.}
#'   \item{schoolf}{Years of schooling of respondent's father.}
#'   \item{schoolm}{Years of schooling of respondent's mother.}
#'   \item{siblings}{Number of siblings.}
#'   \item{wages}{Hourly wages in dollars.}
#'   \item{asvab02}{Arithmetic reasoning.}
#'   \item{asvab03}{Word knowledge.}
#'   \item{lnwages}{Natural logarithm of wages.}
#'   \item{exp2}{Work experience squared.}
#'  }
#' @docType data
#' @usage data(wages2)
#' @keywords datasets
#' @examples
#' data(wages2)
"wages2"

#' Effects of Temperature on Water Consumption
#'
#' An experimental study of the  effects of temperature on water consumption
#'   through self-reported thirst.
#'   Fifty participants were in a room for 4 hours doing a variety of tasks.
#'   Before the experiment,
#'   each participant was acclimated to a standard temperature
#'   of 70 degrees Fahrenheit.
#'
#' @format A dataframe with 50 rows and three variables
#' \describe{
#'   \item{temp}{Room temperature in degrees Fahrenheit.
#'     The temperature was manipulated such that the participants were exposed
#'     to a specific temperature in the room for the 4 hours of the experiment.}
#'   \item{thirst}{Self-reported measure of thirst at the end of a 2-hour period. The scale is from 1 (not at all thirsty) to 5 (very thirsty)}
#'   \item{water}{The number of deciliters of water consumed during the last 2 hours of the study.}
#'  }
#' @docType data
#' @usage data(water)
#' @keywords datasets
#' @examples
#' data(water)
"water"

#' Teacher Expectancies on Student Achievement
#'
#' @format A dataframe with 40 rows and four variabels
#' \describe{
#'   \item{exp}{Teacher expectancy based on an intelligence test given to the student the previous year.}
#'   \item{soc}{The average observer rating of social warmth.}
#'   \item{inp}{The average observer rating of input to the student.}
#'   \item{ach}{Studen achievement measured using a score in the test at the end of the semester.}
#'  }
#' @docType data
#' @usage data(achievement)
#' @keywords datasets
#' @examples
#' data(achievement)
"achievement"

#' Stability of Alienation
#'
#' @format A \eqn{6 \times 6} covariance matrix of the following variables
#' \describe{
#'   \item{anomie67}{Anomie measured in 1967.}
#'   \item{powerless67}{Powerlessness measured in 1967.}
#'   \item{anomie71}{Anomie measured in 1971.}
#'   \item{powerless71}{Powerlessness measured in 1971.}
#'   \item{edu}{Education.}
#'   \item{sei}{Duncan's occupation status index (SEI).}
#'  }
#' @docType data
#' @usage data(alienation)
#' @keywords datasets
#' @examples
#' data(alienation)
"alienation"

#' Simple Mediation Model Using Latent Variables Data
#'
#' Simple mediation model data using latent variables.
#'
#' @docType data
#' @usage data(dat_med_simple_lat)
#' @keywords datasets
#' @examples
#' data(dat_med_simple_lat)
"dat_med_simple_lat"

#' Serial Mediation Model Data
#'
#' Serial mediation model data
#'   from a study of presumed media influence.
#'
#' @docType data
#' @usage data(dat_med_serial2)
#' @keywords datasets
#' @examples
#' data(dat_med_serial2)
"dat_med_serial2"
