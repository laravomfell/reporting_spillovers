# This will run the following combinations

for N_P in "2" "5"
do
    for BW_H in "0.4" "0.6"
    do
        for BW_G  in "1.0" "2.0"
        do
            for BW_TREND in "20" "10"
            do
                for BW_DAILY in "0.041666667" "0.020833333" "0.083333333" # 1, 0.5, 2 hrs
                do
		    for BW_WEEKLY in "0.333333333" "0.166666667" "0.666666667" # 8, 4, 16 hrs
		    do
			Rscript run.R --snapshot --numsmooth ${N_P} --bw_daily ${BW_DAILY} --bw_weekly ${BW_WEEKLY} --bw_trend ${BW_TREND} --bw_g ${BW_G} --bw_h ${BW_H} --follow_trig_prob 0.05 --parents_proportion 0.9 --experimentid "test_clustering" 
		    done
                done
            done
        done
    done
done
