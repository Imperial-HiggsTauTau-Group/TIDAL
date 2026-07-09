for algo in DeepTau2018v2p5 PNet UParT; do 
  for wp in 5 6 7; do
    mkdir Draw/plots/tauID_studies/et_method_1/${algo}_${wp}/Combined
    python Draw/scripts/hadd_cp_datacards.py -i Draw/plots/tauID_studies/et_method_1/${algo}_${wp}/Run3_*/cpdecay/et/*.root -o Draw/plots/tauID_studies/et_method_1/${algo}_${wp}/Combined/added_histo_et.root -c et
  done
done