!
!    Input and output channel information
!
!     OUTCH  : output channel number
!     INCH   : terminal input channel number
!     CINCH  : current input channel
!     INPCH  : channel used by the input command
!     RINCH  : channel used by the read command
!     OOUTCH : file output channel number
!     tOUTCH : trace file output channel number
!     rOUTCH : resid file output channel number
!     qOUTCH : time file output channel number
!     tFLAG  : trace file output flag (value set in main program)
!     rFLAG  : resid file output flag (value set in main program)
!     qFLAG  : time file output flag (value set in main program)
!     terse  : If true causes reduced output to OUTCH and no file output
!-----------------------------------------------------------------------
!
      common /chans/ outch,inch,cinch,rinch,inpch,ooutch,toutch,routch, &
                     qoutch,tflag,rflag,qflag,terse
      integer outch,inch,cinch,rinch,inpch,ooutch,toutch,routch,qoutch
      logical tflag,rflag,qflag,terse

