
"sabre" <-
function(model.formula.uni,
         model.formula.bi=NULL,
         model.formula.tri=NULL,
         case,
         alpha = 0.01,
         approximate = 5,
         max.its = 100,
         arithmetic.type = "fast",
         offset = "",
         convergence = 0.00005,
         correlated = "yes",
         left.end.point = NULL,
         right.end.point = NULL,
         first.family = "binomial",
         second.family = "binomial",
         third.family = "binomial",
         first.link.function = "logit",
         second.link.function = "logit",
         third.link.function = "logit",
         first.mass = 12,
         second.mass = 12,
         third.mass = 12,
         ordered = FALSE,
         first.scale = -10000.0,
         second.scale = -10000.0,
         third.scale = -10000.0,
         first.rho = 0.0,
         second.rho = 0.0,
         third.rho = 0.0,
         first.sigma = 1.0,
         second.sigma = 1.0,
         third.sigma = 1.0,
         tolerance = 0.000001,
         equal.scale = FALSE,
         depend = FALSE,
         only.first.derivatives = FALSE,
	 adaptive.quad = FALSE,
         fixed.effects=FALSE)
{



# some preprocessing

if(fixed.effects && (!is.null(model.formula.bi) || !is.null(model.formula.tri)))
{
   stop("error : cannot have multilevel fixed effects model")
}


if(fixed.effects && (first.family != "gaussian"))
{
   stop("error : cannot have non Gaussian fixed effects model")
}
  
# check endpoints
left.end.point.set<-FALSE
right.end.point.set<-FALSE
# left endpoint
if(!is.null(left.end.point))
{
  # applies only to binary response (binomial) or poisson models
  if(first.family == "guassian" )
   {
     stop("error : left.end.point can only be specified for binomial or Poisson models")
   }
  # check the values
  if(!is.numeric(left.end.point)) 
    {
      stop("error : left.end.point value must be numeric")
    }	
  if(length(left.end.point) > 1)
    {
        stop("error : left.end.point has more than one value")
    }
  if(left.end.point < 0.0)
    {
       stop("error : left.end.point value must be positive")
    } 
   left.end.point.set<-TRUE
}
# right end point
if(!is.null(right.end.point))
{
  # applies only to binary response (binomial) models
  if(first.family != "binomial")
   {
     stop("error : right.end.point  can only be specified for binomial models")
   }
  if(!is.numeric(right.end.point)) 
    {
      stop("error : right.end.point value must be numeric")
    }	
  if(length(right.end.point) > 1)
    {
        stop("error : right.end.point has more than one value")
    }
  if(right.end.point < 0.0)
    {
       stop("error : right.end.point value must be positive")
    }
  right.end.point.set<-TRUE
}
# check the number of variates
if(left.end.point.set || right.end.point.set)
{
  # applies only to univariate binomial models
  if(!is.null(model.formula.bi) || !is.null(model.formula.tri))
   {
     stop("error : endpoints can only be specified for univariate models")
   }
}




if(equal.scale && depend)
{
  stop("cannot have equal.scale == depend == TRUE")
}

if(!is.list(case))
{
   case<-list(case)
}


#         script.file<-"sabreR.sab"
#         data.file<-"sabreR.dat"
#         log.file<-"sabreR.log"

         script.file<-tempfile()
         data.file<-tempfile()
         log.file<-tempfile()

if(ordered)
{
   ordered<-"yes"
}
else
{
  ordered<-"no"
}
if(first.family=="ordered")
{
   ordered<-""
}


# do some sanity checks
if(length(case) == 0)
{
  stop("no case variable specified")
}


# removed this test pre-empting the succesful implementation of AQ for multilevel
# models - DJG 01/02/08
# if(length(case) > 1 && adaptive.quad)
# {
#   stop("adaptive quadrature not availbale for multilevel models")
# }



#if( (length(case) > 1)  && ( !is.null(model.formula.bi) ||  !is.null(model.formula.tri)) && !(depend || equal.scale))
#{
#  stop("multilevel multivariate modelling not supported")
#}


if(length(case) > 2)
{
  stop("maximum number of levels is 3")
}

 # helper functions

sort.data.frame <- function(x, key, ...) {
    if (missing(key)) {
        rn <- rownames(x)
        if (all(rn %in% 1:nrow(x))) rn <- as.numeric(rn)
        x[order(rn, ...), , drop=FALSE]
    } else {
        x[do.call("order", c(x[key], ...)), , drop=FALSE]
    }
}

  
  foldl<-function(binary.op,start.value,values=NULL)
  {
    if(length(values) == 0) # nothing to do
      {
        return(start.value)
      }
    return(foldl(binary.op,binary.op(start.value,values[[1]]),tail(values,length(values)-1)))
  }

  
variable.occurance<-function(model.names)
  {   
    number.of.models<-length(model.names)
    all.names<-foldl(union,model.names[[1]],tail(model.names,number.of.models-1))
    appearances<-vector("list",length(all.names)) # record the occurences of each
    for(name.index in 1:length(all.names))
      {
        occurences<-0
        occurs.in<-NULL
        for(model.index in 1:number.of.models)
          {
            if(all.names[[name.index]] %in% model.names[[model.index]])
              {
                occurences<-occurences+1
                occurs.in<-c(occurs.in,model.index)
              }
          }
        # find out which ones it is in
        
        appearances[[name.index]]<-list(name=all.names[[name.index]],occurences=occurences,occurs.in=occurs.in)
      }
    return(appearances)
  }


sabre.data<-function(models,case)
  {
    # constructs a data frame suitable for use with sabre
    # assumes the models appear in order (ie univariate, then bivariate and so on)
    # number of models
    number.of.models<-length(models)
    # number of observations in each model
    number.of.observations<-NULL  
    model.names<-list()
    model.matrices<-list()
    for(model in models)
      {
        # check the models are consistent
        current.model.frame<-model.frame(model)

        # account for any offset variables
        # get the column numbers of the model frame that contain the offset variable
       offset.vars.col.nums<-attr(attr(current.model.frame,"terms"),"offset")


        number.of.observations<-c(number.of.observations,length( attr(current.model.frame,"row.names")))
        # build the design matrix
        current.model.matrix<-model.matrix(model,current.model.frame)

        # augment the model matrix with any offset variables
        if(!is.null(offset.vars.col.nums))
        {
          for(offset.var.col.num in 1:length(offset.vars.col.nums))
          {
            current.model.matrix<-cbind(current.model.matrix,current.model.frame[offset.vars.col.nums[offset.var.col.num]])
            current.model.matrix<-as.matrix(current.model.matrix)
          }
        }


        # add it to the list
        model.matrices<-c(model.matrices,list(current.model.matrix))
        # and add the names to the list of names
        model.names<-c(model.names,list(attr(current.model.matrix,"dimnames")[[2]]))
      }
    # process the occurences
    occurences<-variable.occurance(model.names)
   # start constructing the data frame
    model.data<-NULL
    if(number.of.models > 1) 
      {
        # need an rvar
        # model.data<-data.frame(r.variate=rep(1:number.of.models,each=number.of.observations))
        model.data<-data.frame(r.variate=rep(1:number.of.models,times=number.of.observations))
      }
   # add the appropriate column for each variate
   for(occurence in occurences)
     {
       new.column<-NULL
       # get a model with this variable in
       model.index<-occurence$occurs.in[[1]]
       for(model.number in occurence$occurs.in) # need a column for each occurence
         {
           new.column.name<-NULL
           if(number.of.models > 1)
             {
               # need to add extra information to the variable names
               new.column.name<-paste(occurence$name,".",sep="")
               new.column.name<-paste(new.column.name,as.character(model.number),sep="")
             }
           else
             {
               new.column.name<-occurence$name
             }
           # create this column - remember it is zero for everything other than the model we are considering
           if(model.number == 1)
             {
               # can just copy it
               new.column<-data.frame(dummy.name=model.matrices[[model.number]][,occurence$name])
             }
           else
             {
               # use zeros
               new.column<-data.frame(dummy.name=rep(0,number.of.observations[1]))
             }
           # change its name
           attr(new.column,"names")[[1]]<-new.column.name
           if(number.of.models > 1)
             {
               for(level in 2:number.of.models)
                 {
                   # now keep adding to it appropriatly
                   new.sub.column<-NULL
                   if(model.number == level)
                     {
                       # can just copy it
                       new.sub.column<-data.frame(dummy.name=model.matrices[[model.number]][,occurence$name])
                     }
                   else
                     {
                       # use zeros
                       new.sub.column<-data.frame(dummy.name=rep(0,number.of.observations[level]))
                     }
                   attr(new.sub.column,"names")[[1]]<-new.column.name
                   # add this level
                   new.column<-rbind(new.column,new.sub.column)
                 } # level
             }
           if(is.null(model.data))
             {
               model.data<-new.column
             }
           else
             {
               model.data<-cbind(model.data,new.column)
             }
         } # occurence
     }
    # tidy up a bit
    rm(model.matrices)
    


    # create the case variables

#if(!(depend || equal.scale))
#{  
#    for(case.var.number in 1:length(case))
#     {
#       case.var<-data.frame(dummy.name=case[case.var.number])
#       if(number.of.models > 1)
#        {
#          for(level in 2:number.of.models)
#            {
#              case.var<-rbind(case.var,case.var)
#            }
#        }
#       case.var.name<-paste("case.",case.var.number,sep="")
#       attr(case.var,"names")[[1]]<-case.var.name
#       cat("first\n")
#       model.data<-cbind(model.data,case.var)
#       cat("end first\n")
#     }
#}
    if(!(depend || equal.scale))
      {
   # is the data unbalanced ?
        balanced.data<-TRUE
        if(length(case) > 1)
          {
            if(length(case[[1]]) != length(case[[2]]))
              {
                balanced.data<-FALSE
              }
          }
        if(balanced.data == TRUE)
          {
            for(case.var.number in 1:length(case))
              {
                case.var<-data.frame(dummy.name=case[case.var.number])
                if(number.of.models > 1)
                  {
                    for(level in 2:number.of.models)
                      {
                        case.var<-rbind(case.var,case.var)
                      }
                  }
                case.var.name<-paste("case.",case.var.number,sep="")
                attr(case.var,"names")[[1]]<-case.var.name
                model.data<-cbind(model.data,case.var)
              }
          }
        else # unalanced data
          {
            model.data<-cbind(model.data,rbind(data.frame(case.1=case[[1]]),data.frame(case.1=case[[2]])))
          }   
      }    
    else
      {
                                        # this is a state dependent model -> does not have multiple levels
        if(length(case) != 2)
          {
            stop("need to specify 2 case variables for state dependent model")
          }
        model.data<-cbind(model.data,rbind(data.frame(case.1=case[[1]]),data.frame(case.1=case[[2]])))
      }
    
    
    # create the response variables
    # the first response
    response.var<-model.frame(models[[1]])[[1]]
    response.vars<-data.frame(dummy.name=response.var)
    if(number.of.models > 1)
      {
        for(level in 2:number.of.models)
          {
            response.var<-model.frame(models[[level]])[[1]]
            response.vars<-rbind(response.vars,data.frame(dummy.name=response.var))
          }
      }
    # tidy up
    rm(response.var)
    # change the name
    attr(response.vars,"names")[[1]]<-"response"
    model.data<-cbind(model.data,response.vars)
    # tidy up
    rm(response.vars)
    
    return(model.data)
  }

    # process the model formula

  # put the models into a list
models<-list(model.formula.uni)
if(!is.null(model.formula.bi))
    {
      models<-c(models,model.formula.bi)
    }
    if(!is.null(model.formula.tri))
    {
      models<-c(models,model.formula.tri)
    }


  # set the number of variates
  model<-NULL
  if(length(models)==1)
    {
      model="univariate"
    }
  else if(length(models)==2)
   {
     model="bivariate"
   }
  else if(length(models)==3) 
    {
      model="trivariate"
    }

   # find out what the offset variables are (if any) - they have to be excluded from the fit command
   offset.variable.names<-NULL
   for(model.formula.number in 1:length(models))
    {
      # create the model frame for the model formula
      current.model.frame<-model.frame(models[[model.formula.number]])
      offset.vars.col.nums<-attr(attr(current.model.frame,"terms"),"offset")
      if(!is.null(offset.vars.col.nums))
      {
        # in the current implementation of sabre can only have one offset variable
	if(length(offset.vars.col.nums) > 1)
          {
              stop("error : current version of sabre only supports a single offset variable")	
          }
	if(model.formula.number > 1)
          {
              stop("error : current version of sabre does not support offset variables in bivariate or trivariate model formula")	
          }
         for(offset.var.col.num in 1:length(offset.vars.col.nums))
          {
    	     offset.variable.names<-c(offset.variable.names,names(current.model.frame)[offset.vars.col.nums[offset.var.col.num]])
          }
      }	
    }

    # cleanse the offset variable names
    if(!is.null(offset.variable.names))
     {
       cleansed.names<-offset.variable.names
       if(length(cleansed.names) > 0)
         {
          for(index in 1:length(cleansed.names))
           {
             cleansed.names[[index]]<-gsub(" ","",cleansed.names[[index]])
           }
        offset.variable.names<-cleansed.names
       }
       rm(cleansed.names)
    }


    # the sabre.data helper function now accounts for offset variables
    model.design.matrix<-sabre.data(models,case)

    if((length(case) == 1) || depend || equal.scale)
      {
        model.design.matrix<-sort.data.frame(model.design.matrix,key="case.1")    
      }
    
    else if(length(case) == 2)
      {


        # is the data unbalanced ?
        balanced.data<-TRUE
        if(length(case[[1]]) != length(case[[2]]))
          {
            balanced.data<-FALSE
          }
      
          if(balanced.data == TRUE)
          {

          model.design.matrix<-sort.data.frame(model.design.matrix,key="case.2")    
          new.model.design.matrix<-data.frame()
          while(nrow(model.design.matrix) != 0)
           {
              sorted.group<-subset(model.design.matrix,model.design.matrix$case.2==model.design.matrix$case.2[1])
              sorted.group<-sort.data.frame(sorted.group,key="case.1")
              model.design.matrix<-subset(model.design.matrix,model.design.matrix$case.2!=model.design.matrix$case.2[1])
              new.model.design.matrix<-rbind(new.model.design.matrix,sorted.group)
           }  
        model.design.matrix<-new.model.design.matrix
        # tidy up
        rm(new.model.design.matrix)
          }
          else
         {
        model.design.matrix<-sort.data.frame(model.design.matrix,key="case.1")    
         }
  
}




  # remove any white space from the variable names    
    cleansed.names<-names(model.design.matrix)
    if(length(cleansed.names) > 0)
      {
        for(index in 1:length(cleansed.names))
         {
           cleansed.names[[index]]<-gsub(" ","",cleansed.names[[index]])
         }
      attr(model.design.matrix,"names")<-cleansed.names
      }
    rm(cleansed.names)
    model.design.matrix.names<-attributes(model.design.matrix)$names
    model.exp.variables<-list()
    for(exp.var.name in model.design.matrix.names)
      {
        if(exp.var.name != "case.1" && exp.var.name != "case.2"  && exp.var.name != "response" && exp.var.name != "r.variate" && !(exp.var.name %in% offset.variable.names))
          {
            # make sure there is no constant in an ordered response model
            if(first.family != "ordered" || 
               (exp.var.name != "(Intercept)" && exp.var.name != "(Intercept).1" &&  exp.var.name != "(Intercept).2"))
            {
                   model.exp.variables<-c(model.exp.variables,exp.var.name)
            }
          }
      }
    model.y.variate<-"response"
  
    # stucuture of model.design.matrix is now [exp.vars] + y.var + case + <r.var>
    # write it to the data file
    write.table(model.design.matrix,file=data.file,col.names=FALSE,row.names=FALSE,quote=FALSE)
    # the names for the data command
    data.command.arguments<-names(model.design.matrix)


    # see if we have constants - if so - set them appropriately
    # defaults
    first.const.var<-""
    second.const.var<-""
    third.const.var<-""
    if("(Intercept)" %in% names(model.design.matrix))
      {
        first.const.var<-"(Intercept)"
      } 
    if("(Intercept).1" %in% names(model.design.matrix))
      {
        first.const.var<-"(Intercept).1"
      } 
    if("(Intercept).2" %in% names(model.design.matrix))
      {
       second.const.var<-"(Intercept).2"
      } 
    if("(Intercept).3" %in% names(model.design.matrix))
      {
       third.const.var<-"(Intercept).3"
      } 

   
    # set the r.variate if required
    # default
    r.variate<-NULL
    if(!is.null(model.formula.bi))
      {
        r.variate<-"r.variate"
      }


  
"process.log" <- function(sabre.model)
  {
    # read the file
    log<-readLines(con=sabre.model$log.file)

    # matches
    iteration.match<-grep("Iteration",log)
    parameter.match<-grep("Parameter",log)
    setting.match<-grep("Setting",log)
    correlation.match<-grep("Correlation",log)
    xvars.match<-grep("X-vars",log)
    name.match<-grep("Name",log)
    # find the prompts <S> 
    prompt.match<-grep("<S>",log)
    
    sabre.model$lfit.estimates.print.message<-paste(log[parameter.match[1]:(prompt.match[prompt.match>parameter.match[1]][1]-1)])
    sabre.model$fit.estimates.print.message<-paste(log[parameter.match[2]:(prompt.match[prompt.match>parameter.match[2]][1]-1)])
    sabre.model$settings.print.message<-paste(log[setting.match[1]:(prompt.match[prompt.match>setting.match[1]][1]-1)])
    try(
        sabre.model$lfit.correlations.print.message<-paste(log[(correlation.match[1]+2):(prompt.match[prompt.match>correlation.match[1]][1]-1)]),
        silent=TRUE)

    try(
        sabre.model$fit.correlations.print.message<-paste(log[(correlation.match[2]+2):(prompt.match[prompt.match>correlation.match[2]][1]-1)]),silent=TRUE
        )
    sabre.model$lfit.iterations.print.message<-paste(log[iteration.match[1]:(prompt.match[prompt.match>iteration.match[1]][1]-1)])
    sabre.model$fit.iterations.print.message<-paste(log[iteration.match[3]:(prompt.match[prompt.match>iteration.match[3]][1]-1)])
    sabre.model$model.print.message.shm<-paste(log[xvars.match[1]:(prompt.match[prompt.match>xvars.match[1]][1]-1)])
    sabre.model$model.print.message.rem<-paste(log[xvars.match[2]:(prompt.match[prompt.match>xvars.match[2]][1]-1)])
    sabre.model$variables.print.message<-paste(log[name.match[1]:(prompt.match[prompt.match>name.match[1]][1]-1)])

    sabre.model
  }
  

"sabre.fit" <- function(sabre.model,
                        # data,
                        script.file="sabreR.sab",
                        data.file="sabreR.dat",
                        log.file="sabreR.log",
                        sabre.binary)
 {

   
#   write.table(subset(data,select=all.variables),file=data.file,col.names=FALSE,row.names=FALSE,quote=FALSE)  

   # where the outout goes
   sabre.script<-"output "
   sabre.script<-c(sabre.script,log.file)
   sabre.script<-c(sabre.script,"\n")
   
   # create a string containing the sabre script
   # the data
   sabre.script<-c(sabre.script,"data ")


   for(data.command.argument in head(data.command.arguments,n=length(data.command.arguments)-1))
     {
       sabre.script<-c(sabre.script,data.command.argument)
       sabre.script<-c(sabre.script," &");
       sabre.script<-c(sabre.script,"\n")
     }
   sabre.script<-c(sabre.script,tail(data.command.arguments,n=1))
   sabre.script<-c(sabre.script,"\n")

   # read the data file
   sabre.script<-c(sabre.script,"read ")
   sabre.script<-c(sabre.script,data.file)
   sabre.script<-c(sabre.script,"\n")
   
   
   # case
   if(length(case) == 1 || sabre.model$depend || sabre.model$equal.scale)
     {
       sabre.script<-c(sabre.script,"case case.1")
       sabre.script<-c(sabre.script,"\n")
     }
   else
     {
       # check for unbalanced data

       balanced.data<-TRUE
       if(length(case) > 1)
         {
           if(length(case[[1]]) != length(case[[2]]))
             {
               balanced.data<-FALSE
             }
         }
       if(balanced.data == TRUE)
         {
           sabre.script<-c(sabre.script,"case first = case.1 second = case.2")
         }
       else
         {
           sabre.script<-c(sabre.script,"case case.1")
           sabre.script<-c(sabre.script,"\n")
         }
       sabre.script<-c(sabre.script,"\n")
     }
   # y.variate
   sabre.script<-c(sabre.script,"yvar ",model.y.variate)
   sabre.script<-c(sabre.script,"\n")
   # r.variate
   if(!is.null(r.variate))
     {
       sabre.script<-c(sabre.script,"rvar ",r.variate)
       sabre.script<-c(sabre.script,"\n")
     }
   # constant variables


   if(sabre.model$first.const.var != "")
     {
       sabre.script<-c(sabre.script,"constant first=",sabre.model$first.const.var)
       sabre.script<-c(sabre.script,"\n")
     }
   if(sabre.model$second.const.var != "")
     {
       sabre.script<-c(sabre.script,"constant second=",sabre.model$second.const.var)
       sabre.script<-c(sabre.script,"\n")
     }
   if(sabre.model$third.const.var != "")
     {
       sabre.script<-c(sabre.script,"constant third=",sabre.model$third.const.var)
       sabre.script<-c(sabre.script,"\n")
     }
   # alpha
   sabre.script<-c(sabre.script,"alpha ",sabre.model$alpha)
   sabre.script<-c(sabre.script,"\n")
   # arithmetic
   sabre.script<-c(sabre.script,"arithmetic ",sabre.model$arithmetic)
   sabre.script<-c(sabre.script,"\n")
   # convergence
   sabre.script<-c(sabre.script,"convergence ",sabre.model$convergence)
   sabre.script<-c(sabre.script,"\n")
   # correlated
   sabre.script<-c(sabre.script,"correlated ",sabre.model$correlated)
   sabre.script<-c(sabre.script,"\n")
   # end points
   if(left.end.point.set && right.end.point.set)
    {
      sabre.script<-c(sabre.script,"endpoints both ",sabre.model$left.end.point," ",sabre.model$right.end.point)
      sabre.script<-c(sabre.script,"\n")
    }
   else if(left.end.point.set)
    {
      sabre.script<-c(sabre.script,"endpoints left ",sabre.model$left.end.point)
      sabre.script<-c(sabre.script,"\n")
    }
   else if(right.end.point.set)
    {
      sabre.script<-c(sabre.script,"endpoints right ",sabre.model$right.end.point)
      sabre.script<-c(sabre.script,"\n")
    }
   else
    {
      sabre.script<-c(sabre.script,"endpoints no")
      sabre.script<-c(sabre.script,"\n")
    }

   # sabre.script<-c(sabre.script,"endpoints ",sabre.model$endpoints)
   # sabre.script<-c(sabre.script,"\n")
   # first family
   sabre.script<-c(sabre.script,"family first= ",sabre.model$first.family)
   sabre.script<-c(sabre.script,"\n")
   # second family
   sabre.script<-c(sabre.script,"family second= ",sabre.model$second.family)
   sabre.script<-c(sabre.script,"\n")
   # third family
   sabre.script<-c(sabre.script,"family third= ",sabre.model$third.family)
   sabre.script<-c(sabre.script,"\n")
   # first link
   sabre.script<-c(sabre.script,"link first= ",sabre.model$first.link.function)
   sabre.script<-c(sabre.script,"\n")
   # second link
   sabre.script<-c(sabre.script,"link second= ",sabre.model$second.link.function)
   sabre.script<-c(sabre.script,"\n")
   # third link
   sabre.script<-c(sabre.script,"link third= ",sabre.model$third.link.function)
   sabre.script<-c(sabre.script,"\n")
   # first mass
   sabre.script<-c(sabre.script,"mass first= ",sabre.model$first.mass)
   sabre.script<-c(sabre.script,"\n")
   # second mass
   sabre.script<-c(sabre.script,"mass second= ",sabre.model$second.mass)
   sabre.script<-c(sabre.script,"\n")
   # third mass
   sabre.script<-c(sabre.script,"mass third= ",sabre.model$third.mass)
   sabre.script<-c(sabre.script,"\n")
   # maximum iterations
   sabre.script<-c(sabre.script,"maxits ",sabre.model$max.its)
   sabre.script<-c(sabre.script,"\n")
   # model
   if(!sabre.model$depend)
    {
      sabre.script<-c(sabre.script,"model ",sabre.model$model)
      sabre.script<-c(sabre.script,"\n")
    }
   # offset
   if(!is.null(offset.variable.names))
   {
      # note - current version of sabre only allows for one offset variable
      sabre.script<-c(sabre.script,"offset ",offset.variable.names[1])
      sabre.script<-c(sabre.script,"\n")
   }	
   # sabre.script<-c(sabre.script,"offset ",sabre.model$offset)
   # sabre.script<-c(sabre.script,"\n")
   # ordering
   sabre.script<-c(sabre.script,"ordered ",sabre.model$ordered)
   sabre.script<-c(sabre.script,"\n")
   # rho
   # sabre.script<-c(sabre.script,"rho ",sabre.model$rho)
   # sabre.script<-c(sabre.script,"\n")
   # rho12
   # sabre.script<-c(sabre.script,"rho 1,2= ",sabre.model$rho12)
   # sabre.script<-c(sabre.script,"\n")
   # rho13
   # sabre.script<-c(sabre.script,"rho 1,3= ",sabre.model$rho13)
   # sabre.script<-c(sabre.script,"\n")
   # rho23
   # sabre.script<-c(sabre.script,"rho 2,3= ",sabre.model$rho23)
   #sabre.script<-c(sabre.script,"\n")

   # first rho
   sabre.script<-c(sabre.script,"rho first=",sabre.model$first.rho)
   sabre.script<-c(sabre.script,"\n")
   # second rho
   sabre.script<-c(sabre.script,"rho second=",sabre.model$second.rho)
   sabre.script<-c(sabre.script,"\n")
   # third rho
   sabre.script<-c(sabre.script,"rho third=",sabre.model$third.rho)
   sabre.script<-c(sabre.script,"\n")


   # first scale
   sabre.script<-c(sabre.script,"scale first= ",sabre.model$first.scale)
   sabre.script<-c(sabre.script,"\n")
   # second scale
   sabre.script<-c(sabre.script,"scale second= ",sabre.model$second.scale)
   sabre.script<-c(sabre.script,"\n")
   # third scale
   sabre.script<-c(sabre.script,"scale third= ",sabre.model$third.scale)
   sabre.script<-c(sabre.script,"\n")
   # first sigma
   sabre.script<-c(sabre.script,"sigma first= ",sabre.model$first.sigma)
   sabre.script<-c(sabre.script,"\n")
   # second sigma
   sabre.script<-c(sabre.script,"sigma second= ",sabre.model$second.sigma)
   sabre.script<-c(sabre.script,"\n")
   # third sigma
   sabre.script<-c(sabre.script,"sigma third= ",sabre.model$third.sigma)
   sabre.script<-c(sabre.script,"\n")
   # tolerance
   sabre.script<-c(sabre.script,"tolerance ",sabre.model$tolerance)
   sabre.script<-c(sabre.script,"\n")


   if(sabre.model$equal.scale)
     {
       sabre.script<-c(sabre.script,"eqscale yes")
       sabre.script<-c(sabre.script,"\n")
     }
   if(sabre.model$depend)
     {
       sabre.script<-c(sabre.script,"depend yes")
       sabre.script<-c(sabre.script,"\n")
     }
   if(sabre.model$only.first.derivatives)
     {
       sabre.script<-c(sabre.script,"der1 yes")
       sabre.script<-c(sabre.script,"\n")
     }
   if(sabre.model$adaptive.quad)
     {
       sabre.script<-c(sabre.script,"quad a")
       sabre.script<-c(sabre.script,"\n")
     }	
   else
    {
       sabre.script<-c(sabre.script,"quad g")
       sabre.script<-c(sabre.script,"\n")
    }	


   # do the fitting
   # the explanatory variables
   #all.exp.variables<-NULL
   all.exp.variables<-model.exp.variables
   #for(exp.variable in model.exp.variables)
   # for(exp.variable in all.variables)
    #  {
        # if(!is.factor(data$exp.variable)) # drop exp vars used as factors
      #   if(!(exp.variable %in% variables.as.factors)) # drop exp vars used as factors
     #     {
            
           # if(exp.variable != sabre.model$case) # drop the case
           #   {
           #     if(exp.variable != sabre.model$y.variate) # drop the y variate
           #       {
      #              all.exp.variables<-c(all.exp.variables,exp.variable) 
           #       }
           #   }
      #    }
     # }
   # all.exp.variables<-c(all.exp.variables,factor.variables) # the names of the factors
   # the constants
   #if(sabre.model$first.const.var != "")
   #  {
   #    all.exp.variables<-c(all.exp.variables,sabre.model$first.const.var) 
   #  }
   #if(sabre.model$second.const.var != "")
   #  {
   #    all.exp.variables<-c(all.exp.variables,sabre.model$second.const.var)
   #  }
   #if(sabre.model$third.const.var != "")
   #  {
   #    all.exp.variables<-c(all.exp.variables,sabre.model$third.const.var)
   #  }
   # lfit and nvars
   # okay - now they have to be in some bloody order !!!!!!
   nvar.first.length<-NULL
   nvar.second.length<-NULL
   if(!is.null(model.formula.bi))
     {
          ordered.all.exp.variables<-NULL
          # find the positions of the first lot
          search.result<-grep("\\.1",all.exp.variables)
          for(index in search.result)
            {
              ordered.all.exp.variables<-c(ordered.all.exp.variables,all.exp.variables[index])
            }   
	   if(length(search.result) > 0)
               {   
	         # first nvars
        	 sabre.script<-c(sabre.script,"nvar first= ",length(search.result))
                 nvar.first.length<-length(search.result)
                 sabre.script<-c(sabre.script,"\n")
               }
          search.result<-grep("\\.2",all.exp.variables)
          for(index in search.result)
            {
              ordered.all.exp.variables<-c(ordered.all.exp.variables,all.exp.variables[index])
            }   
	   if(length(search.result) > 0 && !is.null(model.formula.tri))
               {   
	         # second nvars
        	 sabre.script<-c(sabre.script,"nvar second= ",length(search.result))
                 nvar.second.length<-length(search.result)
                 sabre.script<-c(sabre.script,"\n")
               }   
          search.result<-grep("\\.3",all.exp.variables)
          for(index in search.result)
            {
              ordered.all.exp.variables<-c(ordered.all.exp.variables,all.exp.variables[index])
            }      
	all.exp.variables<-ordered.all.exp.variables
        # tidy up
        rm(ordered.all.exp.variables)	
     }

   sabre.script<-c(sabre.script,"lfit ")
   for(exp.variable in head(all.exp.variables,n=length(all.exp.variables)-1))
     {
       sabre.script<-c(sabre.script,exp.variable)
       sabre.script<-c(sabre.script," &");
       sabre.script<-c(sabre.script,"\n")
     }
   sabre.script<-c(sabre.script,tail(all.exp.variables,n=1))
   sabre.script<-c(sabre.script,"\n")
   sabre.script<-c(sabre.script,"display e")
   sabre.script<-c(sabre.script,"\n")
   sabre.script<-c(sabre.script,"display m")
   sabre.script<-c(sabre.script,"\n")
   # specify - the nvars again !!!!
   if(!is.null(nvar.first.length))
    {
      sabre.script<-c(sabre.script,"nvar first= ",nvar.first.length)
      sabre.script<-c(sabre.script,"\n")
    }
   if(!is.null(nvar.second.length))
    {
      sabre.script<-c(sabre.script,"nvar second= ",nvar.second.length)
      sabre.script<-c(sabre.script,"\n")
    }

   # fit
   # sabre.script<-c(sabre.script,"fit ")
   if(fixed.effects)
    {
      sabre.script<-c(sabre.script,"fefit ")   
    }
   else
   {
     sabre.script<-c(sabre.script,"fit ")   
   }


   for(exp.variable in head(all.exp.variables,n=length(all.exp.variables)-1))
     {
       sabre.script<-c(sabre.script,exp.variable)
       sabre.script<-c(sabre.script," &");
       sabre.script<-c(sabre.script,"\n")
     }
   sabre.script<-c(sabre.script,tail(all.exp.variables,n=1))
   sabre.script<-c(sabre.script,"\n")
   sabre.script<-c(sabre.script,"display e")
   sabre.script<-c(sabre.script,"\n")
   sabre.script<-c(sabre.script,"display m")
   sabre.script<-c(sabre.script,"\n")
   sabre.script<-c(sabre.script,"display s")
   sabre.script<-c(sabre.script,"\n")
   sabre.script<-c(sabre.script,"display v")
   sabre.script<-c(sabre.script,"\n")
   # stop
   sabre.script<-c(sabre.script,"stop")
   sabre.script<-c(sabre.script,"\n")

   
   # write it to file
   cat(as.character(sabre.script),file=script.file,sep="")
   # create the command
   # command.string<-paste(sabre.binary,"<",script.file)<
   #command.string<-paste(sabre.binary," -terse -input=./",script.file,sep="")	
   #system(command.string)
   script.file<-as.character(script.file)
   .Fortran("fortfunc",script.file,nchar(script.file),PACKAGE="sabreR")



 }


  
  # set any conditional defaults
  if(first.scale == -10000.0)
    {
      if(first.family == "poisson")
        {
          first.scale = 0.5
        }
      else
        {
          first.scale = 1.0
        }
    }
    if(second.scale == -10000.0)
    {
      if(second.family == "poisson")
        {
          second.scale = 0.5
        }
      else
        {
          second.scale = 1.0
        }
    }
    if(third.scale == -10000.0)
    {
      if(third.family == "poisson")
        {
          third.scale = 0.5
        }
      else
        {
          third.scale = 1.0
        }
    }



# if necessary - determine the location of the sabre executable from the package structure
sabre.binary<-NULL
if(class(sabre.binary) == "NULL")
  {
    pac.lst<-installed.packages()
    loc<-try(pac.lst["sabreR",2],silent=TRUE)
    if(class(loc) == "try-error")
      {
        stop("cannot find a sabre executable")
      }
    sabre.binary<-paste(loc,"/sabreR/exec/sabre",sep="")
  }


# create the sabre model structure
  sabre.instance<-list(alpha = alpha,
                       approximate = approximate,
                       max.its = max.its,
                       arithmetic.type = arithmetic.type,
                       case = case,
                       offset = offset,
                       convergence = convergence,
                       correlated = correlated,
                       left.end.point = left.end.point,
                       right.end.point = right.end.point,
                       first.family = first.family,
                       second.family = second.family,
                       third.family = third.family,
                       first.link.function = first.link.function,
                       second.link.function = second.link.function,
                       third.link.function = third.link.function,
                       first.mass = first.mass,
                       second.mass = second.mass,
                       third.mass = third.mass,
                       ordered = ordered,
                       first.scale = first.scale,
                       second.scale = second.scale,
                       third.scale = third.scale,
                       first.rho = first.rho,	
                       second.rho = second.rho,
                       third.rho = third.rho,
                       first.sigma = first.sigma,
                       second.sigma = second.sigma,
                       third.sigma = third.sigma,
                       tolerance = tolerance,
                       model= model,
                       first.const.var=first.const.var,
                       second.const.var=second.const.var,
                       third.const.var=third.const.var,
                       script.file=script.file,
                       data.file=data.file,
                       log.file=log.file,
                       equal.scale = equal.scale,
                       depend = depend,
                       only.first.derivatives=only.first.derivatives,
                       adaptive.quad=adaptive.quad,
                       fixed.effects = fixed.effects,
                       model.print.message=NULL,
                       variables.print.message=NULL,
                       lfit.iterations.print.message=NULL,
                       fit.iterations.print.message=NULL,
                       lfit.estimates.print.message=NULL,
                       fit.estimates.print.message=NULL,
                       settings.print.message=NULL,
                       lfit.correlations.print.message=NULL,
                       fit.correlations.print.message=NULL)
  
  attr(sabre.instance,"class")<-"sabre"

  # process the model
  sabre.fit(sabre.instance,sabre.binary=sabre.binary,script.file=script.file,data.file=data.file,log.file=log.file)
  #sabre.fit(sabre.instance,data,sabre.binary=sabre.binary)
  sabre.instance<-process.log(sabre.instance)

   # tidy up
   # command.string<-paste("rm",script.file)
   # system(command.string)
   # command.string<-paste("rm",data.file)
   # system(command.string)
   # command.string<-paste("rm",log.file)
   # system(command.string)
   unlink(script.file)
   unlink(data.file)
   unlink(log.file)

  # return the result
  return(sabre.instance)
}



