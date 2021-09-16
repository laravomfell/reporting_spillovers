@echo off
SetLocal


set NP=2 5 10
set BW_H=0.6 0.8 1.0 1.2
set BW_G=1.0
set BW_TREND=5 10 20
set BW_WEEKLY=0.33333
set BW_DAILY=0.041667


::Rscript run.R --snapshot --numsmooth %%n --bw_h %%h --bw_g %%g --bw_daily %%d --bw_weekly %%w --bw_trend %%t --real_data --experimentid "real_data_main_run" 
::Rscript run.R --snapshot --numsmooth %%n --bw_h %%h --bw_g %%g --bw_daily %%d --bw_weekly %%w --bw_trend %%t --follow_trig_prob 0.05 --parents_proportion 0.9 --experimentid "test_clustering"

for %%n in (%NP%) do ( 
  for %%h in (%BW_H%) do (
    for %%g in (%BW_G%) do (
      for %%t in (%BW_TREND%) do (
        for %%w in (%BW_WEEKLY%) do (
          for %%d in (%BW_DAILY%) do (
    	      echo  %%n %%g %%h %%t %%w %%d
            Rscript run.R --snapshot --numsmooth %%n --bw_h %%h --bw_g %%g --bw_daily %%d --bw_weekly %%w --bw_trend %%t --real_data --experimentid "real_data_main_run"
            Rscript run.R --snapshot --g_delay --numsmooth %%n --bw_h %%h --bw_g %%g --bw_daily %%d --bw_weekly %%w --bw_trend %%t --real_data --experimentid "real_data_main_run"
          )
        )
      )
    )   
  )
)
