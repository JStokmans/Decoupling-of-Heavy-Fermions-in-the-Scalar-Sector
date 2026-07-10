* This FORM code performs a matching scale power counting (powers of M) of the Beetle Diagram in the different regions of integration space.
* This code loads the expressions for the Forests of the Beetle from the same directory as where this code is located.
* These expressions can be generated using the file "Beetle_Diagram_Forests.frm".
* Note that the contribution from the measure has NOT been taken into account and should be corrected for manually.

Dimension 4;
Symbols a,n,m,M;
Symbols [lf^2-M^2],[(lf+p)^2-M^2],
        [(lf+l1)^2-M^2],[(lf+l2)^2-M^2],
        [(lf+p-l1)^2-M^2],[(lf+p-l2)^2-M^2],
        [(lf+p-l1+l2)^2-M^2],[(lf+l1-l2)^2-M^2],
        [(l2-l1)^2-m^2],[(l1-p)^2-m^2],[l2^2-m^2];      * Denominators of I_Gamma
Symbols [(lf-l2)^2-M^2];                                * New denominators of t_gamma_3*I_Gamma
Symbols [(Qf+l1)^2-M^2],[(Qf+l2)^2-M^2],
        [(Qf+Q-l1)^2-M^2],[(Qf+Q-l2)^2-M^2],
        [(Qf+Q-l1+l2)^2-M^2],[(Qf+l1-l2)^2-M^2],
        [(l1-Q)^2-m^2],[(Qf-l2)^2-M^2];                 * New denominators of t_gamma_1*t_gamma_3*t_gamma_f*I_Gamma
Symbols [(lf-l1)^2-M^2],[(lf-l1+l2)^2-M^2],[l1^2-m^2];  * New denominators of t_Gamma*t_gamma_1*t_gamma_3*t_gamma_f*I_Gamma
CFunctions F (Symmetric);
Vectors lf,l1,l2,p;                                     * Unfixed momenta
Vectors Q,Qf;                                           * gamma_1 renormalization scale
.global

* Turning warnings off to prevent FORM from giving one every time a term does not fit inside the bracket function "F" (which happens all the time).
OFF Warnings;
OFF Statistics;

#procedure ScaleFermionDenoms
    id [lf^2-M^2]^n?                = (a^2)^n;
    id [(lf+p)^2-M^2]^n?            = (a^2)^n;
    id [(lf+l1)^2-M^2]^n?           = (a^2)^n;
    id [(lf+l2)^2-M^2]^n?           = (a^2)^n;
    id [(lf+p-l1)^2-M^2]^n?         = (a^2)^n;
    id [(lf+p-l2)^2-M^2]^n?         = (a^2)^n;
    id [(lf+p-l1+l2)^2-M^2]^n?      = (a^2)^n;
    id [(lf+l1-l2)^2-M^2]^n?        = (a^2)^n;
    id [(lf-l2)^2-M^2]^n?           = (a^2)^n;
    id [(Qf+l1)^2-M^2]^n?           = (a^2)^n;
    id [(Qf+l2)^2-M^2]^n?           = (a^2)^n;
    id [(Qf+Q-l1)^2-M^2]^n?         = (a^2)^n;
    id [(Qf+Q-l2)^2-M^2]^n?         = (a^2)^n;
    id [(Qf+Q-l1+l2)^2-M^2]^n?      = (a^2)^n;
    id [(Qf+l1-l2)^2-M^2]^n?        = (a^2)^n;
    id [(Qf-l2)^2-M^2]^n?           = (a^2)^n;
    id [(lf-l1)^2-M^2]^n?           = (a^2)^n;
    id [(lf-l1+l2)^2-M^2]^n?        = (a^2)^n;
#endprocedure

load Beetle_Diagram_Forests.sav;

G [dummy] = [I_G] - [t_g_f] - [t_g_3] + [t_g_3t_g_f] - [t_g_1]
        - [t_G] + [t_Gt_g_f] + [t_Gt_g_3] - [t_Gt_g_3t_g_f] + [t_Gt_g_1];
.sort
hide;

L [R_G] = [dummy];
.sort
hide;

* Power counting around the origin of integration space
G [count 0] = [R_G];
id M = a*M;
id Q = a*Q;
id Qf = a*Qf;
id [(l1-Q)^2-m^2]^n?            = (a^2)^n;
#call ScaleFermionDenoms

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort

Collect F;
id F(a?) = 1;

Print;
.sort
hide [count 0];

* Scaling lf by a factor 'a' to do the power counting in the region {lf}.
G [count lf] = [R_G];
id lf = a*lf;
id M = a*M;
id Q = a*Q;
id Qf = a*Qf;
id [(l1-Q)^2-m^2]^n?            = (a^2)^n;
#call ScaleFermionDenoms

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort

Collect F;
id F(a?) = 1;

Print;
.sort
hide [count lf];


* Scaling l1 by a factor 'a' to do the power counting in the region {l1}.
G [count l1] = [R_G];
id l1 = a*l1;
id M = a*M;
id Q = a*Q;
id Qf = a*Qf;
id [(l1-Q)^2-m^2]^n?            = (a^2)^n;
#call ScaleFermionDenoms

id [(l1-p)^2-m^2]^n?            = (a^2)^n;
id [(l2-l1)^2-m^2]^n?           = (a^2)^n;
id [l1^2-m^2]^n?                = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort

Collect F;
id F(a?) = 1;

Print;
.sort
hide [count l1];


* Scaling l2 by a factor 'a' to do the power counting in the region {l2}.
L [count l2] = [R_G];
id l2 = a*l2;
id M = a*M;
id Q = a*Q;
id Qf = a*Qf;
id [(l1-Q)^2-m^2]^n?            = (a^2)^n;
#call ScaleFermionDenoms

id [(l2-l1)^2-m^2]^n?           = (a^2)^n;
id [l2^2-m^2]^n?                = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort

Collect F;
id F(a?) = 1;

Print;
.sort
hide [count l2];


* Scaling all loop momenta by a factor 'a' to do the power counting in the region {lf,l1}.
G [count lfl1] = [R_G];
id lf = a*lf;
id l1 = a*l1;
id M = a*M;
id Q = a*Q;
id Qf = a*Qf;
id [(l1-Q)^2-m^2]^n?            = (a^2)^n;
#call ScaleFermionDenoms

id [(l1-p)^2-m^2]^n?            = (a^2)^n;
id [(l2-l1)^2-m^2]^n?           = (a^2)^n;
id [l1^2-m^2]^n?                = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort

Collect F;
id F(a?) = 1;

Print;
.sort
hide [count lfl1];



* Scaling all loop momenta by a factor 'a' to do the power counting in the region {lf,l2}.
L [count lfl2] = [R_G];
id lf = a*lf;
id l2 = a*l2;
id M = a*M;
id Q = a*Q;
id Qf = a*Qf;
id [(l1-Q)^2-m^2]^n?            = (a^2)^n;
#call ScaleFermionDenoms

id [(l2-l1)^2-m^2]^n?           = (a^2)^n;
id [l2^2-m^2]^n?                = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort

Collect F;
id F(a?) = 1;

Print;
.sort
hide [count lfl2];



* Scaling all loop momenta by a factor 'a' to do the power counting in the region {l1,l2}.
L [count l1l2] = [R_G];
id l1 = a*l1;
id l2 = a*l2;
id M = a*M;
id Q = a*Q;
id Qf = a*Qf;
id [(l1-Q)^2-m^2]^n?            = (a^2)^n;
#call ScaleFermionDenoms

id [(l1-p)^2-m^2]^n?            = (a^2)^n;
id [(l2-l1)^2-m^2]^n?           = (a^2)^n;
id [l2^2-m^2]^n?                = (a^2)^n;
id [l1^2-m^2]^n?                = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort

Collect F;
id F(a?) = 1;

Print;
.sort
hide [count l1l2];



* Scaling all loop momenta by a factor 'a' to do the power counting in the region {lf,l1,l2}.
L [count lfl1l2] = [R_G];
id lf = a*lf;
id l1 = a*l1;
id l2 = a*l2;
id M = a*M;
id Q = a*Q;
id Qf = a*Qf;
id [(l1-Q)^2-m^2]^n?            = (a^2)^n;
#call ScaleFermionDenoms

id [(l1-p)^2-m^2]^n?            = (a^2)^n;
id [(l2-l1)^2-m^2]^n?           = (a^2)^n;
id [l2^2-m^2]^n?                = (a^2)^n;
id [l1^2-m^2]^n?                = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort

Collect F;
id F(a?) = 1;

Print;
.sort
Delete storage;
.end
