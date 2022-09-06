#!/usr/bin/env nextflow

nextflow.enable.dsl=2
nextflow.preview.recursion=true

// A prototypical active learning workflow for PILs ======================================
//
// This is meant to be a prototype for the active learning workflow in TIPS.
//
// See the parameters sections below for the parameters that can be tuned for
// the workflow, the results will be saved in a `$proj` folder set via params.
//
// TODOs:
// - [ ] modularize & document the subworkflows (train/sample/label/check)
// - [ ] test this with cheaper single point calculations (e.g. lj/lammp)
// - [ ] sanity check for crictical params: [ens_size, geo_size, sp_points]
//
//                                          written by  Yunqi Shao, first ver. 2022.Aug.29
//========================================================================================

// Initial Configraitons =================================================================
params.proj         = 'exp/ekf-v2'
params.restart_from = false
params.init_geo     = 'skel/init/*.xyz'
params.init_model   = 'skel/pinn/pinet-ekf.yml'
params.init_ds      = 'datasets/pils-40ps.{yml,tfr}'
params.init_time    = 0.5
params.init_steps   = 40000
params.ens_size     = 3
params.geo_size     = 6
params.sp_points    = 10
//========================================================================================

// Imports (publish directories are set here) ============================================
include { aseMD } from './tips/nextflow/ase.nf' addParams(publish: "$params.proj/emd")
include { cp2kGenInp } from './tips/nextflow/cp2k.nf' addParams(publish: "$params.proj/cp2k-sp")
include { cp2k } from './tips/nextflow/cp2k.nf' addParams(publish: "$params.proj/cp2k-sp")
include { pinnTrain } from './tips/nextflow/pinn.nf' addParams(publish: "$params.proj/models")
//========================================================================================

// Ietrartion options ====================================================================
params.frmsetol     = 0.150
params.ermsetol     = 0.005
params.fmaxtol      = 2.000
params.emaxtol      = 0.020
params.filters      = "--filter 'peratom(energy)<-125.1' --filter 'abs(force)<15.0'"
params.retrain_step = 20000
params.label_flags  = '-f asetraj --subsample --strategy uniform --nsample 10'
// w. force_std:    = '-f asetraj --subsample --strategy sorted --nsample 10'
params.old_flag     = '--nsample 440'
params.new_flag     = '--psample 100' // 6*10 = 60 pts per iter
params.acc_fac      = 4.0
params.min_time     = 0.5
params.max_time     = 1000.0 // 0.5 2  8 32 128 512 | stop here
params.max_gen      = 20     // 0   1  2  3   4   5
//========================================================================================

// Model specific flags ==================================================================
params.train_flags  = '--log-every 1000 --ckpt-every 10000 --batch 1 --max-ckpts 1 --shuffle-buffer 3000'
params.md_flags     = '--ensemble nvt --T 340 --dt 0.5 --log-every 20'
// params.md_flags     = '--ensemble nvt --T 340 --dt 0.5 --log-every 20 --bias heaviside --kb 1'
params.cp2k_inp     = './skel/cp2k/singlepoint.inp'
params.cp2k_aux     = 'skel/cp2k-aux/*'
//========================================================================================

// Main Iteration and Loops ==============================================================
workflow {
  println("Entering TIPS Active Sampler - Proj. Title: $params.proj")
  // Setup initial training inputs
  init_ds = file(params.init_ds)
  init_geo = file(params.init_geo)
  println("  inital Dataset: ${init_ds.name};")
  println("  ${params.geo_size} geometries to start sampling with;")

  ens_size = params.ens_size.toInteger()
  if (params.restart_from) {
    init_gen = params.restart_from.toString()
    init_models = file("${params.proj}/models/gen${init_gen}/*/model", type:'dir')
    assert ens_size == init_models.size : "ens_size ($ens_size) does not match input ($init_models.size)"
    converge = true
    println("  restarting from gen$init_gen ensemble of size $ens_size;")
    // TODO here one should perhaps remove the old directories
  } else{
    init_gen = '0'
    init_models = file(params.init_model, type:'any')
    if (!(init_models instanceof Path)) {
      assert ens_size == init_models.size() : "ens_size does not match input"
      converge = true
      println("  restarting from an ensemble of size $ens_size;")
    } else {
      println("  starting from scratch with the input $init_models.name of size $ens_size;")
      init_models = [init_models] * ens_size
      converge = false
    }
  }

  // more info about the run
  println("  initial models will ${converge? 'not': ''} be trained before first sampling;")
  println("  sampling for $params.init_time ps, will acc./slow by a factor of $params.acc_fac, minimal $params.min_time ps;")
  println("  training flags: $params.train_flags;")
  println("  sampling flags: $params.md_flags;")
  println("  labelling flags: $params.label_flags;")
  println("  mixing flags for old dataset: $params.old_flag;")
  println("  mixing flags for new dataset: $params.new_flag;")
  println("  models will be extended $params.retrain_step steps per generation;")
  println("  ok, ready to go.")
  steps = params.init_steps.toInteger()
  time = params.init_time.toFloat()
  init_inp = [init_gen, init_geo, init_ds, init_models, steps, time, converge]
  //                 0         1        2            3      4     5         6
  al.recurse(channel.value(init_inp)).until{it[0].toInteger()>params.max_gen||it[5]>params.max_time}
}

// Loop for each iteration =================================================================
workflow al {
  take: ch_inp

  main:
  // retrain or keep the model ============================================================
  ch_inp \
    | view {"[gen${it[0]}] ${it[-1]? 'not training': 'training'} the models."} \
    | branch {gen, geo, ds, models, step, time, converge -> \
              keep: converge
                return [gen, models]
              retrain: !converge
                return [gen, models, ds, (1..params.ens_size).toList(), step]} \
    | set {ch_model}

  ch_model.retrain.transpose(by:[1,3]) \
    | map {gen, model, ds, seed, steps -> \
           ["gen$gen/model$seed", ds, model,
            params.train_flags+
            " --seed $seed --train-steps $steps"+
            (gen.toInteger()==0?' --no-cache':'')]}\
    | pinnTrain

  pinnTrain.out.model \
    | map {name, model -> (name=~/gen(\d+)\/model(\d+)/)[0][1,2]+[model]} \
    | map {gen, seed, model -> [gen, model]} \
    | mix (ch_model.keep.transpose()) \
    | groupTuple(size:params.ens_size) \
    | set {nx_models}
  //=======================================================================================

  // sampling with ensable NN =============================================================
  ch_inp | map {[it[0], it[1], it[5]]} | transpose | set {ch_init_t} // init and time
  nx_models \
    | combine (ch_init_t, by:0)  \
    | map {gen, models, init, t -> \
           ["gen$gen/$init.baseName", models, init, params.md_flags+" --t $t"]} \
    | aseMD
  aseMD.out.traj.set {ch_trajs}
  //=======================================================================================

  // relabel with cp2k ====================================================================
  ch_trajs \
    | map {name, traj -> [name, file(params.cp2k_inp), traj, params.label_flags]} \
    | cp2kGenInp \
    | flatMap {name, inps -> inps.collect {["$name/$it.baseName", it]}} \
    | map {name, inp -> [name, inp, file(params.cp2k_aux)]} \
    | cp2k

  cp2k.out.logs \
    | map {name, logs -> (name=~/(gen\d+\/.+)\/(\d+)/)[0][1,2]+[logs]} \
    | map {name, idx, logs -> [name, idx.toInteger(), logs]} \
    | groupTuple(size:params.sp_points) \
    | mergeDS \
    | set {ch_new_ds}
  //=======================================================================================

  // check convergence ====================================================================
  ch_new_ds \
    | join(ch_trajs) \
    | checkConverge \
    | view {name, geo, msg -> "[$name] ${msg.trim()}"}

  checkConverge.out \
    | map{name,geo,msg-> \
          [(name=~/gen(\d+)\/.+/)[0][1], geo, msg.contains('Converged')]} \
    | groupTuple(size:params.geo_size.toInteger()) \
    | map {gen, geo, conv -> [gen, geo, conv.every()]}
    | set {nx_geo_converge}

  //=======================================================================================

  // mix the new dataset ==================================================================
  ch_inp.map {[it[0], it[2]]}.set{ ch_old_ds }
  ch_new_ds \
    | map {name, idx, ds -> [(name=~/gen(\d+)\/.+/)[0][1], ds]} \
    | groupTuple(size:params.geo_size.toInteger()) \
    | join(ch_old_ds) \
    | map {it+[params.new_flag, params.old_flag]} \
    | mixDS \
    | set {nx_ds}
  //=======================================================================================

  // combine everything for new inputs ====================================================
  ch_inp.map{[it[0], it[4]]}.set {nx_step}
  ch_inp.map{[it[0], it[5]]}.set {nx_time}

  acc_fac = params.acc_fac.toFloat()
  min_time = params.min_time.toFloat()
  retrain_step = params.retrain_step.toInteger()

  nx_geo_converge | join(nx_models) | join(nx_ds) | join(nx_time) | join (nx_step) \
    | map {gen, geo, converge, models, ds, time, step -> \
           [(gen.toInteger()+1).toString(),
            geo, ds, models, \
            converge ? step : step+retrain_step, \
            converge ? time*acc_fac : Math.max(time/acc_fac, min_time), \
            converge]} \
    | view {'-'*80+'\n'+ \
            "[gen${it[0]}] next time scale ${it[5]} ps, ${it[6] ? 'no training planned' : 'next training step '+it[4]}."} \
    | set {nx_inp}
  //=======================================================================================

  emit:
  nx_inp
}

// helper processes ======================================================================
process mixDS {
  tag "gen$gen"
  label 'tips'
  publishDir "$params.proj/mixed/gen$gen"
  input: tuple val(gen), path(newDS, stageAs:'*.traj'), path(oldDS, stageAs:'old/*'), val(newFlag), val(oldFlag)
  output: tuple val(gen), path('mix-ds.{tfr,yml}')

  script:
  """
  tips subsample old/${oldDS[0].baseName}.yml -f pinn -o old-ds -of asetraj $oldFlag
  tips convert ${newDS.join(' ')} -f asetraj -o tmp.traj -of asetraj
  tips subsample tmp.traj -f asetraj -o new-ds -of asetraj $newFlag
  tips convert new-ds.traj old-ds.traj -f asetraj -o mix-ds -of pinn --shuffle $params.filters
  rm {new-ds,old-ds,tmp}.*
  """
}

process mergeDS {
  tag "$name"
  label 'tips'
  publishDir "$params.proj/merge/$name"
  input: tuple val(name), val(idx), path(logs, stageAs:'*.log')
  output: tuple val(name), path('merged.idx'), path('merged.traj')

  script:
  """
  printf "${idx.join('\\n')}" > merged.idx
  tips convert ${logs.join(' ')} -f cp2klog -o merged -of asetraj
  """
}

process checkConverge {
  tag "$name"
  label 'tips'
  publishDir "$params.proj/check/$name"

  input:
  tuple val(name), path(idx), path(logs), path(traj)

  output:
  tuple val(name), path('*.xyz'), stdout

  script:
  fmaxtol = params.fmaxtol
  emaxtol = params.emaxtol
  frmsetol = params.frmsetol
  ermsetol = params.ermsetol
  sp_points = params.sp_points
  """
  #!/usr/bin/env python
  import numpy as np
  from ase import Atoms
  from ase.io import read, write
  from tips.io import load_ds
  from tips.io.filter import filters2fn

  filters = "$params.filters".replace("'", '').split(' ')[1::2]
  filter_fn = filters2fn(filters) # ^ a crude extractor

  idx = [int(i) for i in np.loadtxt("$idx")]
  logs = load_ds("$logs", fmt='asetraj')
  traj = load_ds("$traj", fmt='asetraj')

  idx, logs = tuple(zip(*(
      (i, datum) for (i, datum) in zip(idx, logs) if filter_fn(datum))))

  e_label = np.array([datum['energy']/len(datum['elem']) for datum in logs])
  f_label = np.array([datum['force'] for datum in logs])
  e_pred = np.array([traj[i]['energy']/len(traj[i]['elem']) for i in idx])
  f_pred = np.array([traj[i]['force'] for i in idx])

  ecnt = np.sum(np.abs(e_pred-e_label)>$emaxtol)
  fcnt = np.sum(np.any(np.abs(f_pred-f_label)>$fmaxtol,axis=(1,2)))
  emax = np.max(np.abs(e_pred-e_label))
  fmax = np.max(np.abs(f_pred-f_label))
  ermse = np.sqrt(np.mean((e_pred-e_label)**2))
  frmse = np.sqrt(np.mean((f_pred-f_label)**2))
  converged = (emax<$emaxtol) and (fmax<$fmaxtol) and (ermse<$ermsetol) and (frmse<$frmsetol) and (len(idx)==$sp_points)

  geoname = "$name".split('/')[1]
  if converged:
      msg = f'Converged; will restart from latest frame.'
      new_geo = logs[np.argmax(idx)]
  else:
      msg = f'energy: {ecnt}/{len(idx)} failed, max={emax:.2f} rmse={ermse:.2f}; '\
            f'force: {fcnt}/{len(idx)} failed, max={fmax:.2f} rmse={frmse:.2f}.'
      new_geo = logs[np.argmin(idx)]
  atoms = Atoms(new_geo['elem'], positions=new_geo['coord'], cell=new_geo['cell'],
                pbc=True)
  write(f'{geoname}.xyz', atoms)
  print(msg)
  """
}
//==================================================================== END OF THE WORKFLOW
