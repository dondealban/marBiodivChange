# List of random effects for LME models
rand1 <- ~ 1 | studyName
rand2 <- ~ 1 | studySub
rand3 <- ~ 1 | subSiteID
rand4 <- ~ year0Z | subSiteID

randList1 <- list(rand1, rand2, rand4)
