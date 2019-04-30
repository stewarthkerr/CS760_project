gold_at:
		./generate_gold.py -i ../data/Arab.Meristem/arabidopsis.meristem.expression.csv -c "Feature_ID" --delimiter=',' -o ../data/Arab.Meristem/arabidopsis.meristem.biogrid.interactions.csv --pickle ../data/Arab.Meristem/arabidopsis.meristem.biogrid.graph.pickle

gold_atd:
		python3 -m pdb ./generate_gold.py -i ../data/Arab.Meristem/arabidopsis.meristem.expression.csv -c "Feature_ID" --delimiter=',' -o ../data/Arab.Meristem/arabidopsis.meristem.biogrid.interactions.csv --pickle ../data/Arab.Meristem/arabidopsis.meristem.biogrid.graph.pickle

gold_cpib:
		./generate_gold.py -i ../data/CPIB_Lateral_Root_timecourse/E-MTAB-2565-A-AFFY-2-normalized-expressions.tsv -o ../data/CPIB_Lateral_Root_timecourse/CPIB_Lateral_Root.biogrid.interactions.csv --pickle ../data/CPIB_Lateral_Root_timecourse/CPIB_Lateral_Root_timecourse.biogrid.graph.pickle

gold_cpibd:
		python3 -m pdb ./generate_gold.py -i ../data/CPIB_Lateral_Root_timecourse/E-MTAB-2565-A-AFFY-2-normalized-expressions.tsv -o ../data/CPIB_Lateral_Root_timecourse/CPIB_Lateral_Root.biogrid.interactions.csv --pickle ../data/CPIB_Lateral_Root_timecourse/CPIB_Lateral_Root_timecourse.biogrid.graph.pickle

cluster2csv:
	 R CMD BATCH cluster2csv.R

golden2cluster_at:
		./golden2cluster.py -i ../data/Arab.Meristem/arabidopsis.meristem.expression.csv -c "Feature_ID" --delimiter=',' --cluster ../output/Clustering.csv -o ../data/Arab.Meristem/arabidopsis.meristem.modules.interactions.csv --pickle ../data/Arab.Meristem/arabidopsis.meristem.biogrid.graph.pickle --opickle ../data/Arab.Meristem/arabidopsis.meristem.modules.graph.pickle

golden2cluster_atd:
		python3 -m pdb ./golden2cluster.py -i ../data/Arab.Meristem/arabidopsis.meristem.expression.csv -c "Feature_ID" --delimiter=',' --cluster ../output/Clustering.csv -o ../data/Arab.Meristem/arabidopsis.meristem.modules.interactions.csv --pickle ../data/Arab.Meristem/arabidopsis.meristem.biogrid.graph.pickle --opickle ../data/Arab.Meristem/arabidopsis.meristem.modules.graph.pickle
