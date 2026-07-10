* This code computes the expressions for each term in the Forest Formula for the Rocket Diagram and puts these over a common denominator so they can be added together later.
* It stores the resulting expressions in a file "Rocket_Diagram_Forests.sav", which is placed in the same directory as this code.
* Note that if the computation takes an unreasonably long time, one can choose to comment out one of the distinct permutations of the Rocket Diagram. 

Symbols m,M;
Symbols a,n;
Symbols [l1^2-M^2],[(l1-l)^2-M^2],[(l1+p1+p2)^2-M^2],[(l1+p2)^2-M^2],
        [(l1+l)^2-M^2],[(l1+p1+l)^2-M^2],[(l1-p2)^2-M^2],
        [(l+p1+p2)^2-m^2],[l^2-m^2],
        [l2^2-M^2],[(l2+p1+p2)^2-M^2],[(l2-l)^2-M^2];       * Denominators of the Rocket
Vectors p1,p2,l,l1,l2;      * unfixed momenta
Functions f (Symmetric);
CFunctions g (Symmetric);
Vectors x,y,z;
Indices i,j;

OFF Statistics;

*** Definitions of procedures we will need later.
#procedure FactorDenominators
Bracket 
    [l1^2-M^2],
    [(l1-l)^2-M^2],
    [(l1+p1+p2)^2-M^2],
    [(l1+p2)^2-M^2],
    [(l1+l)^2-M^2],
    [(l1+p1+l)^2-M^2],
    [(l1-p2)^2-M^2],
    [(l+p1+p2)^2-m^2],
    [l^2-m^2],
    [l2^2-M^2],
    [(l2+p1+p2)^2-M^2],
    [(l2-l)^2-M^2];
#endprocedure

#define step "1"
#define MaxStep "12"
* Defining a procedure that writes the current progress to the terminal.
* This can then be downloaded to check the progress when running this code on a cluster.
* The procedure requires a preprocessor variable 'step' to be defined beforehand.
* To use it: call the procedure BEFORE each noteworthy (computationally expensive) step.
#procedure ProgressReport
    #message Step 'step'/'MaxStep'
    #redefine step "{'step'+1}"
#endprocedure

*** Defining the proper forests
* Atomic forests
L [I_Gamma] =   (
                2*
                ((g_(1,l1)+M*gi_(1))*(g_(1,l1)-g_(1,l)+M*gi_(1))*(g_(1,l1)+g_(1,p1)+g_(1,p2)+M*gi_(1))*(g_(1,l1)+g_(1,p2)+M*gi_(1))) / 
                ([l1^2-M^2]*[(l1-l)^2-M^2]*[(l1+p1+p2)^2-M^2]*[(l1+p2)^2-M^2]) 

                +((g_(1,l1)+M*gi_(1))*(g_(1,l1)-g_(1,p2)+M*gi_(1))*(g_(1,l1)+g_(1,p1)+g_(1,l)+M*gi_(1))*(g_(1,l1)+g_(1,l)+M*gi_(1))) / 
                ([l1^2-M^2]*[(l1-p2)^2-M^2]*[(l1+p1+l)^2-M^2]*[(l1+l)^2-M^2])
                ) * 

                1/([(l+p1+p2)^2-m^2]*[l^2-m^2]) * 

                ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+g_(2,p1)+g_(2,p2)+M*gi_(2))*(g_(2,l2)-g_(2,l)+M*gi_(2))) / 
                ([l2^2-M^2]*[(l2+p1+p2)^2-M^2]*[(l2-l)^2-M^2]);

L [t_gamma_1*I_Gamma] = (
                        3*
                        ((g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))) / 
                        ([l1^2-M^2]^4)
                        ) *

                        1/([(l+p1+p2)^2-m^2]*[l^2-m^2]) * 

                        ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+g_(2,p1)+g_(2,p2)+M*gi_(2))*(g_(2,l2)-g_(2,l)+M*gi_(2))) / 
                        ([l2^2-M^2]*[(l2+p1+p2)^2-M^2]*[(l2-l)^2-M^2]);

L [t_gamma_2*I_Gamma] = (
                        2*
                        ((g_(1,l1)+M*gi_(1))*(g_(1,l1)-g_(1,l)+M*gi_(1))*(g_(1,l1)+g_(1,p1)+g_(1,p2)+M*gi_(1))*(g_(1,l1)+g_(1,p2)+M*gi_(1))) / 
                        ([l1^2-M^2]*[(l1-l)^2-M^2]*[(l1+p1+p2)^2-M^2]*[(l1+p2)^2-M^2]) 

                        +((g_(1,l1)+M*gi_(1))*(g_(1,l1)-g_(1,p2)+M*gi_(1))*(g_(1,l1)+g_(1,p1)+g_(1,l)+M*gi_(1))*(g_(1,l1)+g_(1,l)+M*gi_(1))) / 
                        ([l1^2-M^2]*[(l1-p2)^2-M^2]*[(l1+p1+l)^2-M^2]*[(l1+l)^2-M^2])
                        ) *

                        1/([(l+p1+p2)^2-m^2]*[l^2-m^2]) * 

                        ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))) / 
                        ([l2^2-M^2]^3);

* Composite forests
L [t_gamma_1*t_gamma_2*I_Gamma] =   (
                                    3*
                                    ((g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))) / 
                                    ([l1^2-M^2]^4)
                                    ) *

                                    1/([(l+p1+p2)^2-m^2]*[l^2-m^2]) * 

                                    ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))) / 
                                    ([l2^2-M^2]^3);

* Non-proper forests
L [t_Gamma*I_Gamma] = (
                2*
                ((g_(1,l1)+M*gi_(1))*(g_(1,l1)-g_(1,l)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))) / 
                ([l1^2-M^2]*[(l1-l)^2-M^2]*[l1^2-M^2]*[l1^2-M^2]) 

                +((g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+g_(1,l)+M*gi_(1))*(g_(1,l1)+g_(1,l)+M*gi_(1))) / 
                ([l1^2-M^2]*[l1^2-M^2]*[(l1+l)^2-M^2]*[(l1+l)^2-M^2])
                ) * 

                1/([l^2-m^2]*[l^2-m^2]) * 

                ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))*(g_(2,l2)-g_(2,l)+M*gi_(2))) / 
                ([l2^2-M^2]*[l2^2-M^2]*[(l2-l)^2-M^2]);

L [t_Gamma*t_gamma_1*I_Gamma] = (
                        3*
                        ((g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))) / 
                        ([l1^2-M^2]^4)
                        ) *

                        1/([l^2-m^2]*[l^2-m^2]) * 

                        ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))*(g_(2,l2)-g_(2,l)+M*gi_(2))) / 
                        ([l2^2-M^2]*[l2^2-M^2]*[(l2-l)^2-M^2]);

L [t_Gamma*t_gamma_2*I_Gamma] = (
                        2*
                        ((g_(1,l1)+M*gi_(1))*(g_(1,l1)-g_(1,l)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))) / 
                        ([l1^2-M^2]*[(l1-l)^2-M^2]*[l1^2-M^2]*[l1^2-M^2]) 

                        +((g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+g_(1,l)+M*gi_(1))*(g_(1,l1)+g_(1,l)+M*gi_(1))) / 
                        ([l1^2-M^2]*[l1^2-M^2]*[(l1+l)^2-M^2]*[(l1+l)^2-M^2])
                        ) *

                        1/([l^2-m^2]*[l^2-m^2]) * 

                        ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))) / 
                        ([l2^2-M^2]^3);

L [t_Gamma*t_gamma_1*t_gamma_2*I_Gamma] = (
                                    3*
                                    ((g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))) / 
                                    ([l1^2-M^2]^4)
                                    ) *

                                    1/([l^2-m^2]*[l^2-m^2]) * 

                                    ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))) / 
                                    ([l2^2-M^2]^3);

* Common denominator
L [Den(R_Gamma)] =
[l2^2-M^2]*
[(l2+p1+p2)^2-M^2]*
[(l2-l)^2-M^2]*
[(l+p1+p2)^2-m^2]*
[l^2-m^2]*
[l1^2-M^2]*

[(l1-l)^2-M^2]*
[(l1+p1+p2)^2-M^2]*
[(l1+p2)^2-M^2]*
[(l1+l)^2-M^2]*
[(l1+p1+l)^2-M^2]*
[(l1-p2)^2-M^2]*
1
;

Trace4,1;
Trace4,2;
.sort
hide;

* Defining the forests over a common denominator using abbreviated notation to facilitate storage.
G [I_G]             = [I_Gamma]/[Den(R_Gamma)];
G [t_g_1]           = [t_gamma_1*I_Gamma]/[Den(R_Gamma)];
G [t_g_2]           = [t_gamma_2*I_Gamma]/[Den(R_Gamma)];
G [t_g_1t_g_2]      = [t_gamma_1*t_gamma_2*I_Gamma]/[Den(R_Gamma)];

G [t_G]             = [t_Gamma*I_Gamma]/[Den(R_Gamma)];
G [t_Gt_g_1]        = [t_Gamma*t_gamma_1*I_Gamma]/[Den(R_Gamma)];
G [t_Gt_g_2]        = [t_Gamma*t_gamma_2*I_Gamma]/[Den(R_Gamma)];
G [t_Gt_g_1t_g_2]   = [t_Gamma*t_gamma_1*t_gamma_2*I_Gamma]/[Den(R_Gamma)];

* Putting things over a common denominator
id [l2^2-M^2]^n?neg_           = [l2^2-M^2]^(4+n)          * g([l2^2-M^2]^-3);
id [(l2+p1+p2)^2-M^2]^n?neg_   = [(l2+p1+p2)^2-M^2]^(2+n)  * g([(l2+p1+p2)^2-M^2]^-1);
id [(l2-l)^2-M^2]^n?neg_       = [(l2-l)^2-M^2]^(2+n)      * g([(l2-l)^2-M^2]^-1);
id [(l+p1+p2)^2-m^2]^n?neg_    = [(l+p1+p2)^2-m^2]^(2+n)   * g([(l+p1+p2)^2-m^2]^-1);
id [l^2-m^2]^n?neg_            = [l^2-m^2]^(3+n)           * g([l^2-m^2]^-2);
id [l1^2-M^2]^n?neg_           = [l1^2-M^2]^(5+n)          * g([l1^2-M^2]^-4);

id [(l1-l)^2-M^2]^n?neg_       = [(l1-l)^2-M^2]^(2+n)      * g([(l1-l)^2-M^2]^-1);
id [(l1+p1+p2)^2-M^2]^n?neg_   = [(l1+p1+p2)^2-M^2]^(2+n)  * g([(l1+p1+p2)^2-M^2]^-1);
id [(l1+p2)^2-M^2]^n?neg_      = [(l1+p2)^2-M^2]^(2+n)     * g([(l1+p2)^2-M^2]^-1);
id [(l1+l)^2-M^2]^n?neg_       = [(l1+l)^2-M^2]^(3+n)      * g([(l1+l)^2-M^2]^-2);
id [(l1+p1+l)^2-M^2]^n?neg_    = [(l1+p1+l)^2-M^2]^(2+n)   * g([(l1+p1+l)^2-M^2]^-1);
id [(l1-p2)^2-M^2]^n?neg_      = [(l1-p2)^2-M^2]^(2+n)     * g([(l1-p2)^2-M^2]^-1);

* Expanding the numerator
#call ProgressReport
id [l1^2-M^2]           = (l1.l1-M^2);
.sort
#call ProgressReport
id [(l1-l)^2-M^2]       = (l1.l1-2*l1.l+l.l-M^2);
.sort
#call ProgressReport
id [(l1+p1+p2)^2-M^2]   = (l1.l1+2*l1.p1+2*l1.p2+p1.p1+2*p1.p2+p2.p2-M^2);
.sort
#call ProgressReport
id [(l1+p2)^2-M^2]      = (l1.l1+2*l1.p2+p2.p2-M^2);
.sort
#call ProgressReport
id [(l1+l)^2-M^2]       = (l1.l1+l.l+2*l1.l-M^2);
.sort
#call ProgressReport
id [(l1+p1+l)^2-M^2]    = (l1.l1+p1.p1+l.l+2*l1.p1+2*l1.l+2*p1.l-M^2);
.sort
#call ProgressReport
id [(l1-p2)^2-M^2]      = (l1.l1+p2.p2-2*l1.p2-M^2);
.sort
#call ProgressReport
id [(l+p1+p2)^2-m^2]    = (l.l+2*l.p1+2*l.p2+p1.p1+2*p1.p2+p2.p2-m^2);
.sort
#call ProgressReport
id [l^2-m^2]            = (l.l-m^2);
.sort
#call ProgressReport
id [l2^2-M^2]           = (l2.l2-M^2);
.sort
#call ProgressReport
id [(l2+p1+p2)^2-M^2]   = (l2.l2+2*l2.p1+2*l2.p2+p1.p1+2*p1.p2+p2.p2-M^2);
.sort
#call ProgressReport
id [(l2-l)^2-M^2]       = (l2.l2-2*l2.l+l.l-M^2);
.sort
#message Expanding numerators complete!

* Restoring the denominators (removing the noncommuting function 'g').
id g(a?) = a;
#message Denominators restored!

#call FactorDenominators

.store
#message Saving results...
save Rocket_Diagram_Forests.sav;
Delete storage;
.end
