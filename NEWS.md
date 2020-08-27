
# v6_9
# 22/06/20 drop ... from argument list to force criterion to pick up information from model object

# 12/09/19 retainModels = FALSE only returns the summary table (useful when models are large)

# 11/09/19 catch for errors in model fits

# 10/09/19 fitModels = FALSE returns the candidate models without fitting them

# 10/09/19 allows for 2D smoothers - see notes for explanation of how to use

# 29/03/17 allows more general likelihood calculation

# 14/02/17 fix bug when smoothTerms is missing; also deal with random effects in the formula (lme4)

# 01/01/17 essential argument allows terms to always be present

# 01/01/17 by argument allows the by terms to be expanded when the by variable is continuous

# 22/12/16 correct bug when picking up AIC for gamm4 models


# 21/12/16 adjust for gamm models


# mygamm4 news

# 11/08/19 add call as attribute of model object
# 10/09/19 mer_only argument to return only the mer object (saves time and space!)
# 11/09/19 allow calc.derivs to be passed in control argument (saves time)

# 1_4
# 16/07/20 ensure start, when supplied, is only given to the first stage of the optimisation (glmer)
