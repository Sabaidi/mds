require(kernlab)
library(kernlab)


binding_data <- read.csv(file="C:\\Users\\micha\\Documents\\GoogleDrive\\5.Semester\\Medical\\Assignments\\03\\BindingData.csv")


peptide_kernel <- function(x,y){
  res = 0
  x_split <- strsplit(x["peptide.sequence"], "")
  y_split <- strsplit(y["peptide.sequence"], "")
  x_split <- unlist(x_split)
  y_split <- unlist(y_split)
  for (i in 1:length(x_split)){
    print(x_split[i])
    if (x_split[i] == y_split[i]){
      res <- res + 1
    }
  }
  return(res)
}

class(peptide_kernel) <- "kernel"


dirac_kernel <- function(x,y){
  if (x["allele"] == y["allele"]){
    return(1)
  }
  else
    return(0)
}

class(dirac_kernel) <- "kernel"


uniform_kernel <- function(x,y){
  return(1)
}

class(uniform_kernel) <- "kernel"

multitask_kernel <- function(x,y){
  return(uniform_kernel(x,y) + dirac_kernel(x,y))
}

class(multitask_kernel) <- "kernel"

combination_kernel <- function(x,y){
  print("Debug")
  return (peptide_kernel(x,y) * multitask_kernel(x,y))
}

class(combination_kernel) <- "kernel"

mat = as.matrix(binding_data)

k = kernelMatrix(kernel = combination_kernel,mat)

Cs = c(10^-4, 10^-3, 10^-2, 10^-1,10^0, 10^1,10^2 ,10^3 ,10^4 )
model = ksvm(binding.label ~ ., data=binding_data, type="C-svc",C=Cs, cross=10, kernel=combination_kernel)
model
