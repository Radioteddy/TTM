  Õ  9   k820309    ?          19.1        ©ü^                                                                                                          
       C:\Users\AkhmetovF\work\Damage project\TTM\ttm_repo\Transfer_Matrix.f90 TRANSFER_MATRIX                                                     
                                                           
                                                           
                         @                                '0                    #POL    #IN_ANGLE    #LAMBDA    #TPULSE    #TPEAK 	   #FLUENCE 
                                                                                                                                   
                                                             
                                                             
                                              	                
                                              
     (          
                     @               @                'h                   #SYSCOMPOUND    #MATERIALS    #LAYER    #SPACEGRID    #ABSORP    #NLAYER    #PHASE    #THETA    #REF    #TRANS    #M    #MLAYER    #FB    #TT    #RR    #T    #R    #A    #POINTS                                                    d                      .                                                     h                             &                                                                                                       °                 
            &                                                                                                ø                 
            &                                                                                               @                
            &                                                                                                                           &                                                                                                Ð                            &                                                                                                                           &                                                                                                `             	               &                                                                                                ¨             
               &                                                                                                ð                            &                   &                                                                                                P                            &                   &                   &                                                                                               È                            &                   &                                                                                           (                                                            8                                                            H         
                                                   P         
                                                   X         
                                                    `                                                             
                 
          #       -DTû!	@        3.1415926535897932384626433832795D0#         @                                                       #SP !   #EP "   #NUMPOINTS #   #ARR $             
                                 !     
                
                                 "     
                
                                  #                                                     $                   
 	              &                                           #         @                                  %                    #ARRAY &   #SUBARRAY '           
                               &                   
               &                                                     
                                '                   
              &                                           #         @                                  (                    #THETA_IN )   #N_IN *   #N_OUT +   #LASER ,   #THETA_OUT -   #R .   #T /             
  @                              )                     
                                 *                     
                                 +                     
                                  ,     0              #SOURCE              D @                              -                      D                                .                      D                                /            #         @                                  0                    #LASER 1   #INTARGET 2             
  @                               1     0              #SOURCE              
D @                               2     h              #MULTILAYER    #         @                                  3                    #LASER 4   #INTARGET 5             
  @                               4     0              #SOURCE              
D @                               5     h              #MULTILAYER    #         @                                   6                    #LASER 7   #INTARGET 8             
  @                               7     0              #SOURCE              
D @                               8     h              #MULTILAYER           `      fn#fn       @   J   CONSTANTS    @  @   J   OBJECTS      @   J   ADD_PROCEDURES    À         SOURCE+OBJECTS #   W  P   a   SOURCE%POL+OBJECTS (   §  H   a   SOURCE%IN_ANGLE+OBJECTS &   ï  H   a   SOURCE%LAMBDA+OBJECTS &   7  H   a   SOURCE%TPULSE+OBJECTS %     H   a   SOURCE%TPEAK+OBJECTS '   Ç  H   a   SOURCE%FLUENCE+OBJECTS #           MULTILAYER+OBJECTS /   '  P   a   MULTILAYER%SYSCOMPOUND+OBJECTS -   w     a   MULTILAYER%MATERIALS+OBJECTS )        a   MULTILAYER%LAYER+OBJECTS -   §     a   MULTILAYER%SPACEGRID+OBJECTS *   ;     a   MULTILAYER%ABSORP+OBJECTS *   Ï     a   MULTILAYER%NLAYER+OBJECTS )   c     a   MULTILAYER%PHASE+OBJECTS )   ÷     a   MULTILAYER%THETA+OBJECTS '   	     a   MULTILAYER%REF+OBJECTS )   
     a   MULTILAYER%TRANS+OBJECTS %   ³
  ¬   a   MULTILAYER%M+OBJECTS *   _  Ä   a   MULTILAYER%MLAYER+OBJECTS &   #  ¬   a   MULTILAYER%FB+OBJECTS &   Ï  H   a   MULTILAYER%TT+OBJECTS &     H   a   MULTILAYER%RR+OBJECTS %   _  H   a   MULTILAYER%T+OBJECTS %   §  H   a   MULTILAYER%R+OBJECTS %   ï  H   a   MULTILAYER%A+OBJECTS *   7  H   a   MULTILAYER%POINTS+OBJECTS             G_PI+CONSTANTS (     p       LINSPACE+ADD_PROCEDURES +     @   a   LINSPACE%SP+ADD_PROCEDURES +   Â  @   a   LINSPACE%EP+ADD_PROCEDURES 2     @   a   LINSPACE%NUMPOINTS+ADD_PROCEDURES ,   B     a   LINSPACE%ARR+ADD_PROCEDURES +   Î  a       CONCATENATE+ADD_PROCEDURES 1   /     a   CONCATENATE%ARRAY+ADD_PROCEDURES 4   »     a   CONCATENATE%SUBARRAY+ADD_PROCEDURES    G         FRESNEL !   Ú  @   a   FRESNEL%THETA_IN      @   a   FRESNEL%N_IN    Z  @   a   FRESNEL%N_OUT      T   a   FRESNEL%LASER "   î  @   a   FRESNEL%THETA_OUT    .  @   a   FRESNEL%R    n  @   a   FRESNEL%T    ®  a       TRANSFERMATRIX %     T   a   TRANSFERMATRIX%LASER (   c  X   a   TRANSFERMATRIX%INTARGET    »  a       AMPL      T   a   AMPL%LASER    p  X   a   AMPL%INTARGET %   È  a       CALCULATE_ABSORPTION +   )  T   a   CALCULATE_ABSORPTION%LASER .   }  X   a   CALCULATE_ABSORPTION%INTARGET 