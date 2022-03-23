#!/bin/bash
# specify paths
log_dir = "."
input_dir=""
output_dir="../results/HLCA_ext_label_transfer/per_dataset"

for i in Kaminski_2020 Eils_2020 Regev_2021_Nuclei Zhang_2021 Meyer_2021_5prime Budinger_2021 Barbry_unpubl Sheppard_2020 Guo_2020_LAM1_3 Sun_2020_batch1 KULeuven_Thienpont_2018Lambrechts_v2 Banovich_Kropski_2020 Sun_2020_batch4 Janssen_2020 Duong_lungMAP_unpubl Lambrechts_2021 Meyer_2021_3prime Misharin_Budinger_2018 Peer_Massague_2020 Gomperts2021_UCLA Wunderink_2021_cryo Schiller_2021 Gomperts_2021_CSMC Lafyatis_Rojas_2019_disease Shalek_2018 Sims_2019 Schiller_2020 Gomperts_2021_CFF Sun_2020_batch3 Wunderink_2021_fresh Regev_2021_Cryo KULeuven_Thienpont_2018Lambrechts_v1 Sun_2020_batch2 Regev_2021_Fresh
do
srun -o ${log_dir}/label_transfer_$i.o.log -e ${log_dir}/label_transfer_$i.e.log -p cpu_p -c 24 -t 24:00:00 --mem=200G --export=ALL conda run -n nn python ./label_transfer.py $i {output_dir} &
done
echo "finished knn classifier for every dataset"
