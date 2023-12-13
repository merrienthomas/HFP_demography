#### 15- EFFECT OF HFP ON DEMOGRAPHIC TRAITS WITH CORRECTIONS ####

selection_with_interaction <- c(2, 14, 3, 15)
selection_without_interaction <- c(2, 14, 3, 15)

# Model for animals with mean HFP value -------------------------------------------------------------
set.seed(15)
output_animal <- Running_MCMCglmm(dfanimal, list_traits_animal, value_checking_animal, list_transform_animal, animal_pop_tree)

# Model for plants with mean HFP value --------------------------------------------------------------
set.seed(16)
output_plant <- Running_MCMCglmm(dfplant, list_traits_plant, value_checking_plant, list_transform_plant, plant_pop_tree)


#### 16-1 - OUTPUTS OF THE SIMULATION ####

write.csv(output_animal[[1]], file = "Animal_phylo_results.csv")
write.csv(output_animal[[2]], file = "Animal_results_correction.csv")
write.csv(output_animal[[3]], file = "Animal_residual_mixed_model.csv")

write.csv(output_plant[[1]], file = "Plant_phylo_results.csv")
write.csv(output_plant[[2]], file = "Plant_results_correction.csv")
write.csv(output_plant[[3]], file = "Plant_residual_mixed_model.csv")


#### 16-2- SIMULATION FOR ANIMALS DEPENDING ON THE MIGRATION DISTANCE ####

for (i in 1:length(dfanimal$SpeciesAccepted)){
  dfanimal$migration[i] <- animal_migration$Distance[animal_migration$SpeciesAccepted==dfanimal$SpeciesAccepted[i]]
}

set.seed(15)
output_animal_long_dist <- Running_MCMCglmm(dfanimal[dfanimal$migration=="Long",], list_traits_animal, value_checking_animal, list_transform_animal, animal_pop_tree)
output_animal_short_dist <- Running_MCMCglmm(dfanimal[dfanimal$migration=="Short",], list_traits_animal, value_checking_animal, list_transform_animal, animal_pop_tree)

write.csv(output_animal_long_dist[[1]], file = "Animal_long_phylo_results.csv")
write.csv(output_animal_long_dist[[2]], file = "Animal_long_results_correction.csv")
write.csv(output_animal_long_dist[[3]], file = "Animal_long_residual_mixed_model.csv")

write.csv(output_animal_short_dist[[1]], file = "Animal_short_phylo_results.csv")
write.csv(output_animal_short_dist[[2]], file = "Animal_short_results_correction.csv")
write.csv(output_animal_short_dist[[3]], file = "Animal_short_residual_mixed_model.csv")

#### 16-3- TRAIT ANALYSIS WITHOUT CORRECTION ####

#Besoin de stocker les résultats des différentes régression dans un seul objet : glmm_stock. 
glmm <- as.data.frame()

for (i in 1:length(list_traits_animal)){
  
  if (value_checking_animal[k] == 1){
    
    list <- Removing_outliers_and_false_values(dfanimal, list_traits_animal[k])
    
    dfbis <- list[[1]]
    storage <- list[[2]]
    
  } else {
    
    if (value_checking[k] == 2){
      
      list <- Removing_outliers(dfanimal, list_traits_animal[k])
      
      dfbis <- list[[1]]
      storage <- list[[2]]
      
    } else {
      
      #if (value_checking[k] == 3){
      
      list <- Removing_outliers_and_zero(dfanimal, list_traits_animal[k])
      
      dfbis <- list[[1]]
      storage <- list[[2]]
      
      #}
      
    }
  }
  
  #How to transform the data
  if (list_transform_animal[k] == 1){
    
    dfbis$transform <- dfbis[, list_traits_animal[k]]
    
  } else {
    
    if (list_transform_animal[k] == 2) {
      
      dfbis$transform <- log(dfbis[, list_traits_animal[k]] + 0.000001)
      
    } else {
      
      if (list_transform_animal[k] == 3) {
        
        dfbis$transform <- sqrt(dfbis[, list_traits_animal[k]] + 0.000001)
        
      } else {
        
        if (list_transform_animal[k] == 4) {
          
          dfbis$transform <- (1/(dfbis[, list_traits_animal[k]] + 0.000001))
          
        } else {
          
          #if (list_transform[k] == 5) {
          
          dfbis$transform <- (sqrt(dfbis[, list_traits_animal[k]] + 0.0001))
          
          #}
          
        }
        
      }
      
    }
    
  }
  
  glmm[k] <- glmer(dfbis[,list_traits_animal[k]] ~ dfbis$coord1dim * dfbis$coord2dim + (1|dfbis$SpeciesAccepted), family = gaussian)
}



