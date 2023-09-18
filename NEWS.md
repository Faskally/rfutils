# rfultils 0.0.3 (16/07/20)

o ensure start, when supplied, is only given to the first stage of the optimisation (glmer)

# rfutils 0.0.2 (22/06/20)

o 22/06/20 drop ... from argument list to force criterion to pick up information from model object
o 12/09/19 retainModels = FALSE only returns the summary table (useful when models are large)
o 11/09/19 catch for errors in model fits
o 10/09/19 fitModels = FALSE returns the candidate models without fitting them
o 10/09/19 allows for 2D smoothers - see notes for explanation of how to use
o 29/03/17 allows more general likelihood calculation
o 14/02/17 fix bug when smoothTerms is missing; also deal with random effects in the formula (lme4)
o 01/01/17 essential argument allows terms to always be present
o 01/01/17 by argument allows the by terms to be expanded when the by variable is continuous
o 22/12/16 correct bug when picking up AIC for gamm4 models
o 21/12/16 adjust for gamm models

# rfutils 0.0.1

o 11/08/19 add call as attribute of model object
o 10/09/19 mer_only argument to return only the mer object (saves time and space!)
o 11/09/19 allow calc.derivs to be passed in control argument (saves time)
