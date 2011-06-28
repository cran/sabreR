

#print.sabre<-function(x,...)
#{

#print.sabre.impl<-function(x,
		           #estimates = TRUE,
                           #iterations = FALSE,
                           #variables = FALSE,
                           #settings = TRUE,
                           #shm = TRUE,
                           #rem = TRUE)

print.sabre<-function(x,
		      estimates = TRUE,
                      iterations = FALSE,
                      variables = FALSE,
                      settings = TRUE,
                      shm = TRUE,
                      rem = TRUE,
                      ...)

  {
    
    # the different print methods
    print.sabre.shm.estimates<-function()
      {
        header<-c("\n","(Standard Homogenous Model)","\n")
        cat(header,sep="\n")
        cat(x$lfit.estimates.print.message,sep="\n")
      }

    print.sabre.rem.estimates<-function()
      {
        header<-c("\n","(Random Effects Model)","\n")
        cat(header,sep="\n")
        cat(x$fit.estimates.print.message,sep="\n")
      }


    print.sabre.shm.iterations<-function()
      {
        header<-c("\n","(Standard Homogenous Model)","\n")
        cat(header,sep="\n")
        cat(x$lfit.iterations.print.message,sep="\n")
      }

    print.sabre.rem.iterations<-function()
      {
        header<-c("\n","(Random Effects Model)","\n")
        cat(header,sep="\n")
        cat(x$fit.iterations.print.message,sep="\n")
      }

    print.sabre.variables<-function()
      {
        cat(x$variables.print.message,sep="\n")
      }


    print.sabre.rem.settings<-function()
      {
        cat(x$model.print.message.rem,sep="\n")
      }

    print.sabre.shm.settings<-function()
      {
        cat(x$model.print.message.shm,sep="\n")
      }


    # the print logic
    
    if(estimates == TRUE)
      {
        if(shm == TRUE)
          {
            print.sabre.shm.estimates()
          }
        if(rem == TRUE)
          {
            print.sabre.rem.estimates()        
          }
      }

    if(settings == TRUE)
      {
        if(shm == TRUE)
          {
            print.sabre.shm.settings()
          }
        if(rem == TRUE)
          {
            print.sabre.rem.settings()        
          }

      #  print.sabre.settings()        
      }
    

    if(iterations == TRUE)
      {
        if(shm == TRUE)
          {
            print.sabre.shm.iterations()
          }
        if(rem == TRUE)
          {
            print.sabre.rem.iterations()        
          }
      }
    
    if(variables == TRUE)
      {
        print.sabre.variables()        
      }

  }
#print.sabre.impl(x,...)
#}