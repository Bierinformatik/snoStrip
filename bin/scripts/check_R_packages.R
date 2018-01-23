pkgTest <- function(x)
  {
    if ( x %in% rownames(installed.packages() ) )
      {
        cat("package found\n")
      }
    else{
      cat( "package not found\n" )
    }
  }




## COMMAND LINE OPTIONS
args <- commandArgs(TRUE)

pkgTest( args[1] )
