for algo in DeepTau2018v2p5 UParT PNet; do
  for wp in 5 6 7; do
    sed -i "s|^output_path:.*|output_path: /vols/cms/dmw25/TIDAL/Draw/plots/tauID_studies/et_method_1/${algo}_${wp}/|" Draw/scripts/cpdecay_datacards_LateRun3.yaml
    sed -i "s|^tau_id:.*|tau_id: ${algo},${wp}|" Draw/scripts/cpdecay_datacards_LateRun3.yaml
    python Draw/scripts/makeDatacards.py --config Draw/scripts/cpdecay_datacards_LateRun3.yaml --batch
  done
done