for algo in DeepTau2018v2p5 PNet UParT; do 
  for wp in 5 6 7; do
    python Draw/scripts/hadd_cp_datacards.py -i Draw/plots/tauID_studies/${algo}_${wp}/Run3_*/cpdecay/et/*.root -o Draw/Plots/TauIDStudyDatacards/${algo}_${wp}/added_histo_et.root -c et
  done
done