####
getWeeklyIncidence <- function(x0, sampleWeeklyIncidence, no_doses, epFlag, A, ep_t, dayNoAfterBurn, weekNo) {
  
  if (dayNoAfterBurn == 0) {
    for (a in 1:A) {
      # ensure r indexing
      baseIndex <- 455 * (a - 1) + 72 * 6
      # Resetting specific elements of x0 to 0
      x0[baseIndex + 2] <- 0
      x0[baseIndex + 3] <- 0
      x0[baseIndex + 4] <- 0
      x0[baseIndex + 5] <- 0
      
      # Resetting a range of elements in x0 to 0
      for (j in 1:9) {
        x0[baseIndex + 8 + j] <- 0
      }
    }
  }
  
  if (dayNoAfterBurn %% 7 == 0 && dayNoAfterBurn > 0) {
    for (a in 1:A) {
      baseIndex <- 455 * (a - 1) + 72 * 6
      for (j in 1:9) {
        incidenceIndex <- baseIndex + 8 + j
        if (epFlag) {
          sampleWeeklyIncidence[weekNo, 9 * (a - 1) + j] <- x0[incidenceIndex] * ep_t[a]
        } else {
          sampleWeeklyIncidence[weekNo, 9 * (a - 1) + j] <- x0[incidenceIndex]
        }
        x0[incidenceIndex] <- 0
      }
      
      doseIndices <- (baseIndex + 2):(baseIndex + 5)
      no_doses[weekNo, 1:4] <- no_doses[weekNo, 1:4] + x0[doseIndices]
      x0[doseIndices] <- 0
    }
    
    weekNo <- weekNo + 1
  }
  
  dayNoAfterBurn <- dayNoAfterBurn + 1
  
  # Return the updated values
  list(x0 = x0, sampleWeeklyIncidence = sampleWeeklyIncidence, no_doses = no_doses, dayNoAfterBurn = dayNoAfterBurn, weekNo = weekNo)
}
