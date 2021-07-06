# These will be run with trend bandwidth 20.

Rscript run.R --regenerate --numsmooth 5   --follow_trig_prob 0.0 --parents_proportion 1.0 
Rscript run.R --regenerate --numsmooth 10  --follow_trig_prob 0.0 --parents_proportion 1.0 
Rscript run.R --regenerate --numsmooth 20  --follow_trig_prob 0.0 --parents_proportion 1.0 
Rscript run.R --regenerate --numsmooth 50  --follow_trig_prob 0.0 --parents_proportion 1.0 
Rscript run.R --regenerate --numsmooth 100 --follow_trig_prob 0.0 --parents_proportion 1.0 


Rscript run.R --regenerate --numsmooth 5   --follow_trig_prob 0.05 --parents_proportion 0.8 
Rscript run.R --regenerate --numsmooth 10  --follow_trig_prob 0.05 --parents_proportion 0.8 
Rscript run.R --regenerate --numsmooth 20  --follow_trig_prob 0.05 --parents_proportion 0.8 
Rscript run.R --regenerate --numsmooth 50  --follow_trig_prob 0.05 --parents_proportion 0.8 
Rscript run.R --regenerate --numsmooth 100 --follow_trig_prob 0.05 --parents_proportion 0.8 
