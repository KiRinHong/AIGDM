#!/bin/bash
find results_omni_final/ -type f -name "Test3*T5.txt" -exec cat {} + > rslts_type1_T5_omni.txt
find results_omni_final/ -type f -name "Test3*T20.txt" -exec cat {} + > rslts_type1_T20_omni.txt
find results_omni_final/ -type f \( -name "Test1*T5.txt" -o -name "Test2*T5.txt" \) -exec cat {} + > rslts_power_T5_omni.txt
find results_omni_final/ -type f \( -name "Test1*T20.txt" -o -name "Test2*T20.txt" \) -exec cat {} + > rslts_power_T20_omni.txt
# find results_omni_new/ -type f \( -name "Test1*T5.txt" -o -name "Test2*T5.txt" \) -exec cat {} + > rslts_power_T5_perturbZero.txt
