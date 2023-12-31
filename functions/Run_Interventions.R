####
# Initialize vectors
populationPerAgeGroup <- numeric()
eta <- numeric()
modelIncidencePerTime <- numeric()
pA <- numeric()
ep_t <- numeric()

# Initialize numeric variables
dailyBirthRate <- 0
totPopulation <- 0
valueLogLikelihood <- 0
currentODETime <- 0
run_start <- 0
run_burn <- 0
run_full <- 0
dt <- 0

# Initialize ageStratification as a numeric vector
ageStratification <- numeric()

# Initialize integer variables
dayNoAfterBurn <- 0
weekNo <- 0
monthNo <- 0

RunInterventions <- function(dailyBirthRate_t, totPopulation_t, ageStratification_t) {
  # variables 
  dailyBirthRate <- dailyBirthRate_t
  totPopulation <- totPopulation_t
  ageStratification <- ageStratification_t
  A <- length(ageStratification)
  eta <- numeric(A)
  populationPerAgeGroup <- numeric(A)
  modelIncidencePerTime <- numeric(A)
  eta[1] <- 0
  
  # Using for loop for calculations
  for (i in 1:(A - 1)) {
    populationPerAgeGroup[i] <- dailyBirthRate * 365 * (ageStratification[i + 1] - ageStratification[i])
    eta[i + 1] <- 1.0 / (365.0 * (ageStratification[i + 1] - ageStratification[i]))
  }
  modelIncidencePerTime[1:(A - 1)] <- 0
  populationPerAgeGroup[A] <- totPopulation - (dailyBirthRate * 365) * ageStratification[A]
  eta[A] <- dailyBirthRate / (totPopulation - (dailyBirthRate * 365) * ageStratification[A])
  modelIncidencePerTime[A] <- 0
  
  # Other variables 
  dt <- 1
  currentODETime <- 0
  dayNoAfterBurn <- 0
  valueLogLikelihood <- 0
  
  # new values
  list(
    A = A,
    eta = eta,
    populationPerAgeGroup = populationPerAgeGroup,
    modelIncidencePerTime = modelIncidencePerTime,
    dt = dt,
    currentODETime = currentODETime,
    dayNoAfterBurn = dayNoAfterBurn,
    valueLogLikelihood = valueLogLikelihood
  )
}

# initialise pVHR, pHR, pLR
pVHR <- numeric() 
pHR <- numeric()   
pLR <- numeric()

