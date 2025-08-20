# PRPP_SMART_BHM
Bayesian hierarchical modeling methods to estimate DTRs in a Partially Randomized, Patient Preference, Sequential, Multiple Assignment, Randomized Trial with a continous end-of-trial outcome. Companion to the manuscript "Bayesian Dynamic Borrowing Approaches for Incorporating Patient Treatment Preferences in SMART Designs" by Sarah Medley, Satrajit Roychoudhury, Thomas M. Braun, and Kelley M. Kidwell. 

## The PRPP-SMART Design
Our PRPP-SMART design assumes there are two stages with two treatment options per stage. At the beginning of stage 1, participants are asked if they have a preference between the two stage 1 treament options (A, B). All participants with a preference are assigned to their preferred treatment while all others are randomized to one of the two treatment options. At the end of stage 1, response status is determined (i.e., responder or non-responder). Responders continue their stage 1 treatment in stage 2 while non-responders are re-assigned to two new treatment options that may be more effective (C, D). Treatment preference is elicited from non-responders at the beginning of stage 2, and non-responders with a preference receive their preferred treatment while all other non-responders are randomized. 

![](Images/PRPP_SMART.png)

PRPP-SMARTs are an alternative to standard SMART designs when strong patient preferences may hinder trial recruitment or retention. PRPP-SMARTs can be visualized as a collection of 4 two-stage sub-trials, one of which is a standard SMART design (panel a).

![](Images/design.png)

Like a standard SMART design, there are 6 *treatment sequences* a participant in the PRPP-SMART may receive: A, AC, AD, B, BC, and BD. A and B correspond to treatment sequences for stage 1 responders while AC, AD, BC, and BC correspond to treatment sequences for stage 1 non-responders. When we incorporate treatment preference into the sequence, there are 20 total *trial pathways* which are denoted in the gray columns of the figure above (notation explained below). 

## PRPP-SMART Notation and Goal
Notation:
- A<sub>1</sub> = {A, B} and A<sub>2</sub> = {C, D} denote stage 1 and 2 treatment, respectively
- T<sub>1</sub> and T<sub>2</sub> are binary indicators of stage 1 and 2 treatment, respectively with T<sub>1</sub> = 1 for A and 0 for B, T<sub>2</sub> = 1 for C and 0 for D
- P<sub>1</sub> and P<sub>2</sub> are binary indicators of stage 1 and 2 preference, respectively

We denote trial pathways by A<sub>1P<sub>1</sub></sub> for responders and A<sub>1P<sub>1</sub></sub>A<sub>2P<sub>2</sub></sub> for non-responders. We observe participant data at the trial pathway level in a PRPP-SMART, but our interest is in estimating dynamic treatment regimens (DTRs) which average over responder and non-responder pathways. We denote DTRs by [A<sub>1</sub>A<sub>1</sub>A<sub>2</sub>]<sub>P<sub>1</sub>P<sub>2</sub></sub>, and there are 16 DTRs subject to treatment preference embedded in our PRPP-SMART design (i.e., the 4 traditional DTRs AAC, AAD, BBC, and BBD but within each preference combination). In particular, we aim to estimate the four indifference DTRs [AAC]<sub>00</sub>, [AAD]<sub>00</sub>, [BBC]<sub>00</sub>, and [BBD]<sub>00</sub> with minimal bias while utilizing all data from a PRPP-SMART. Since the use of non-randomized data makes bias a primary concern, we evaluate our method in a variety of scenarios. 
