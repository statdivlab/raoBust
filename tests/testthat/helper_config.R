# data from `test-offsets`
my_dff <- data.frame("x0" = 1:50,
                     "offs" = 1/(1:50))
my_dff$y0 <- pmax(0.1, rnorm(50, mean = 10*my_dff$x0 * my_dff$offs, sd = 1))

# data from `test-test-statistic`
age <- c(72, 81, 90, 72, 70, 72, 75, 75, 67, 70, 86, 71, 75, 70, 79,
         70, 81, 72, 68, 67, 84, 71, 70, 71, 89, 78, 71, 86, 71, 70, 71,
         70, 73, 73, 67, 80, 74, 75, 70, 72, 81, 70, 71, 70, 69, 82, 74,
         74, 72, 68, 73, 84, 79, 76, 86, 71, 67, 75, 72, 75, 78, 70, 74,
         68, 73, 70, 74, 86, 70, 72, 72, 76, 68, 70, 73, 79, 72, 68, 71,
         67, 71, 68, 81, 83, 73, 79, 69, 82, 77, 84, 81, 83, 80, 76, 75,
         79, 75, 72, 75, 79)
sex <- c("Male", "Female", "Male", "Male", "Female", "Male", "Male",
         "Male", "Female", "Female", "Female", "Male", "Female", "Female",
         "Male", "Male", "Male", "Male", "Male", "Female", "Male", "Male",
         "Male", "Female", "Female", "Male", "Male", "Male", "Male", "Male",
         "Male", "Female", "Male", "Female", "Male", "Female", "Female",
         "Female", "Female", "Female", "Female", "Female", "Female", "Male",
         "Male", "Male", "Female", "Male", "Male", "Male", "Female", "Male",
         "Male", "Male", "Female", "Male", "Male", "Female", "Female",
         "Male", "Female", "Female", "Male", "Male", "Male", "Female",
         "Female", "Male", "Male", "Male", "Female", "Female", "Female",
         "Male", "Female", "Male", "Female", "Male", "Male", "Female",
         "Male", "Female", "Male", "Female", "Male", "Male", "Male", "Male",
         "Female", "Male", "Female", "Male", "Female", "Female", "Male",
         "Female", "Male", "Female", "Male", "Female")
weight <- c(173, 139, 145, 190, 153, 154.5, 161.5, 158.2, 168, 127, 139,
            181.2, 137, 137, 181, 180.5, 151.5, 171, 210, 160, 154, 191,
            181, 124, 111, 172.5, 167, 167, 152.5, 183.5, 204, 135, 208,
            168.5, 230, 88, 163, 86, 167, 167, 158, 192, 126, 155, 239, 146,
            122, 162, 204.2, 168, 136, 156, 170, 193.5, 182, 155, 196, 145.5,
            135, 159, 118, 141, 241, 158.5, 152, 107, 171, 156, 177, 152,
            166, 159, 137, 131, 145, 174, 170, 135, 170.5, 142, 186, 217,
            161, 98, 124, 176, 198, 180.5, 152, 128, 137, 153.5, 130, 97,
            163, 179, 150, 157, 164.5, 139)
mri_mini <- data.frame(age, sex, weight)

nsim <- 200
nn <- 100
simulate_random_data <- function() {
  covariate1 <- rnorm(nn)
  yy <- rpois(nn, lambda = 5 + 0*covariate1)
  df <- data.frame(covariate1 = covariate1,
                   yy = yy)
  return(df)
}

# data from `test-use-case`
my_df <- structure(list(yy = c(87928L, 139438L, 163794L, 53806L, 63075L,
                               27126L, 290005L, 104522L, 170096L, 131047L, 35917L, 36306L, 139066L,
                               247089L, 106406L, 397543L, 35169L, 153151L, 20801L, 36922L),
                        xstar = c(17, 12, 13, 21, 18, 28, 15, 15, 15, 16, 21, 19,
                                  17, 14, 22, 15, 26, 18, 26, 20),
                        xx = c(94, 78, 86, 101,
                               92, 111, 108, 88, 97, 97, 93, 86, 103, 99, 120, 115, 110,
                               110, 99, 90),
                        predictor = c(-1.71008143821379, -1.87180217690159,
                                      -1.88939793879197, -1.57059807911784, -1.63141681915288,
                                      -1.37732569113713, -1.97408102602201, -1.769286613376, -1.86666077740117,
                                      -1.8021222562636, -1.48807705542983, -1.50990831708707, -1.80151564417342,
                                      -1.95606252051933, -1.69644928942373, -2.03688192726104,
                                      -1.44238382777093, -1.81010860789625, -1.33702331211311,
                                      -1.50407739677627)),
                   class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -20L))