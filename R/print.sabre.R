
print.sabre<-function(x,
                      estimates = TRUE,
                      # correlations = FALSE,
                      iterations = FALSE,
                      variables = FALSE,
                      settings = TRUE,
                      shm = TRUE,
                      rem = TRUE,
                      ...)
  {

    sabre.instance<-x    

    # the different print methods
    print.sabre.shm.estimates<-function()
      {
        header<-c("\n","(Standard Homogenous Model)","\n")
        cat(header,sep="\n")
        cat(sabre.instance$lfit.estimates.print.message,sep="\n")
      }

    print.sabre.rem.estimates<-function()
      {

        # header<-c("\n","(Random Effects Model)","\n")
        if(sabre.instance$fixed.effects)
         {
           header<-c("\n","(Fixed Effects Model)","\n")
         }
        else
         {
           header<-c("\n","(Random Effects Model)","\n")           
         }
        cat(header,sep="\n")
        cat(sabre.instance$fit.estimates.print.message,sep="\n")
      }


    #print.sabre.shm.correlations<-function()
    #  {
    #    if(is.null(sabre.instance$lfit.correlations.print.message))
    #      {
    #        cat("correlations are not available for this model\n")
    #      }
    #    else
    #      {
    #        header<-c("\n","Correlation Matrix (Standard Homogenous Model)","_____________________________________________")
    #        cat(header,sep="\n")
    #        cat(sabre.instance$lfit.correlations.print.message,sep="\n")
    #      }
    #  }

    #print.sabre.rem.correlations<-function()
    #  {
    #    if(is.null(sabre.instance$fit.correlations.print.message))
    #      {
    #        cat("correlations are not available for this model\n")
    #      }
    #    else
    #      {
    #        header<-c("\n","Correlation Matrix (Random Effects Model)","________________________________________")
    #        cat(header,sep="\n")
    #        cat(sabre.instance$fit.correlations.print.message,sep="\n")
    #      }
    #  }


    print.sabre.shm.iterations<-function()
      {
        header<-c("\n","(Standard Homogenous Model)","\n")
        cat(header,sep="\n")
        cat(sabre.instance$lfit.iterations.print.message,sep="\n")
      }

    print.sabre.rem.iterations<-function()
      {
        # header<-c("\n","(Random Effects Model)","\n")
        if(sabre.instance$fixed.effects)
         {
           header<-c("\n","(Fixed Effects Model)","\n")
         }
        else
         {
           header<-c("\n","(Random Effects Model)","\n")           
         }
        cat(header,sep="\n")
        cat(sabre.instance$fit.iterations.print.message,sep="\n")
      }

    print.sabre.variables<-function()
      {
        cat(sabre.instance$variables.print.message,sep="\n")
      }

    print.sabre.settings.shm<-function()
      {
        header<-c("\n","(Standard Homogenous Model)","\n")
        cat(header,sep="\n")
        cat(sabre.instance$model.print.message.shm,sep="\n")
      }

    print.sabre.settings.rem<-function()
      {
        # header<-c("\n","(Random Effects Model)","\n")
        if(sabre.instance$fixed.effects)
         {
           header<-c("\n","(Fixed Effects Model)","\n")
         }
        else
         {
           header<-c("\n","(Random Effects Model)","\n")           
         }
        cat(header,sep="\n")
        cat(sabre.instance$model.print.message.rem,sep="\n")
      }


    # the print logic

    if(settings == TRUE)
      {
        if(shm == TRUE)
          {
            print.sabre.settings.shm()
          }
      }
    
    if(estimates == TRUE)
      {
        if(shm == TRUE)
          {
            print.sabre.shm.estimates()
          }
      }


    if(settings == TRUE)
      {
        if(rem == TRUE)
          {
            print.sabre.settings.rem()
          }
      }

    
    if(estimates == TRUE)
      {
        if(rem == TRUE)
          {
            print.sabre.rem.estimates()
          }
      }

    
    # if(correlations == TRUE)
    #  {
    #    if(shm == TRUE)
    #      {
    #        print.sabre.shm.correlations()
    #      }
    #    if(rem == TRUE)
    #      {
    #        print.sabre.rem.correlations()        
    #      }
    #  }

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
