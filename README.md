# PRPP_SMART_BHM
Bayesian hierarchical modeling methods to estimate DTRs in a Partially Randomized, Patient Preference, Sequential, Multiple Assignment, Randomized Trial with a continous end-of-trial outcome. Companion to the manuscript "Bayesian Dynamic Borrowing Approaches for Incorporating Patient Treatment Preferences in SMART Designs" by Sarah Medley, Satrajit Roychoudhury, Thomas M. Braun, and Kelley M. Kidwell. 

## The PRPP-SMART Design
Our PRPP-SMART design assumes there are two stages with two treatment options per stage. At the beginning of stage 1, participants are asked if they have a preference between the two stage 1 treament options (A, B). All participants with a preference are assigned to their preferred treatment while all others are randomized to one of the two treatment options. At the end of stage 1, response status is determined (i.e., responder or non-responder). Responders continue their stage 1 treatment in stage 2 while non-responders are re-assigned to two new treatment options that may be more effective (C, D). Treatment preference is elicited from non-responders at the beginning of stage 2, and non-responders with a preference receive their preferred treatment while all other non-responders are randomized. 

![](Images/PRPP_SMART.png)

