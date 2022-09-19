#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
 This is the workflow for an trial run, it has the following entry points
 `-entry lmpinit` generates the initial geometry guess with moltemplate and lammps
 `-entry cp2k` generates the cp2k trajectories or extends the trajectory
 `-entry train` trains the inital MD models and perform some basic analyses
*/

//======================================//
// THE GENERATION OF INITIAL GEOMETRIES //
//======================================//
include { moltemplate; tag2inp } from './nextflow/build.nf' addParams(publish: 'trajs/build')
include { lammpsMD } from './tips/nextflow/lammps.nf' addParams(publish: 'trajs/lmp')

params.tags = 'a32b32i0,a0b0i32,a16b16i16'
params.rhos = '1.1551,1.0753' // g cm^-3
params.lmp_inp = 'skel/lmp/init.in'

workflow lmpinit {
  channel.fromList(params.tags.tokenize(',')) \
    | combine (channel.fromList(params.rhos.toString().tokenize(','))) \
    | map { tag, rho -> tag2inp(tag, rho) } \
    | moltemplate \
    | map {name, aux -> [name, file(params.lmp_inp), aux]} \
    | lammpsMD
}


//==========================================//
// THE GENERATION OF THE INITAL TRANING SET //
//==========================================//
params.cp2k_inp = 'skel/cp2k/nvt-10ps.inp'
params.cp2k_aux = 'skel/cp2k-aux/*'
params.cp2k_geo = 'trajs/lmp/*/equi.dump'
params.cp2k_from = '10'
include { cp2kMD; cp2k } from './tips/nextflow/cp2k.nf' addParams(publish: 'trajs/cp2k')

workflow cp2kinit {
  channel.fromPath(params.cp2k_geo) \
    | map {it -> \
           ["nvt-0-10ps/${it.parent.name}", file(params.cp2k_inp), it, \
            "--fmt lammps-dump --idx -1 --emap ${file("trajs/build/${it.parent.name}/system.data")}"]} \
    | cp2kMD
}

workflow cp2krestart {
  int from = params.cp2k_from.toInteger()
  channel.fromPath("trajs/cp2k/nvt-${from-10}-${from}ps/*/cp2k-md-1_${2000*from}.restart") \
    | map {it -> \
           ["nvt-${from}-${from+10}ps/${it.parent.name}", it, file(params.cp2k_aux)] }\
    | cp2k
}


//===============================================//
// GENERATION OF INITIAL GEOMETRY FOR PRODUCTION //
//===============================================//

params.lmp2init = 'trajs/lmp/*/equi.dump'

workflow mklmpgeo {
  channel.fromPath(params.lmp2init) \
    | map {it -> \
           ["${it.parent.name}", it, \
            file("trajs/build/${it.parent.name}/system.data")]} \
    | lmp2xyz
}

process lmp2xyz {
  publishDir './skel/init/lmp-geo', mode: 'copy'

  input:
  tuple val(name), path(dump), path(emap)

  output:
  tuple val(name), path('*.xyz')

  script:
  """
  #!/usr/bin/env python
  import re
  from tips.io import load_ds
  m = re.search(r'(a\\d+b\\d+i\\d+)-rho(1\\.\\d+)', '$name')
  shortname =  f'{m[1]}-r{float(m[2]):.2f}'
  ds = load_ds('$dump', fmt='lammps-dump')[-1:].map_elems('$emap')
  ds.convert(f'{shortname}.xyz', fmt='extxyz')
  """
}

params.cp2k2init = 'trajs/cp2k/nvt-40-50ps/*'

workflow mkcp2kgeo {
  channel.fromPath(params.cp2k2init, type:'dir') \
    | map {it -> ["${it.name}", it]} \
    | cp2k2xyz
}

process cp2k2xyz {
  publishDir './skel/init/cp2k-geo', mode: 'copy'

  input:
  tuple val(name), path(dump)

  output:
  tuple val(name), path('*.xyz')

  script:
  """
  #!/usr/bin/env python
  import re
  from tips.io import load_ds
  m = re.search(r'a(\\d+)b(\\d+)+i(\\d+)-rho(1\\.\\d+)', '$name')
  assert m[1]==m[2]
  shortname =  f'm{m[1]}i{m[3]}-r{float(m[4]):.2f}'
  ds = load_ds('$dump/cp2k-md', fmt='cp2k')[-1:]
  ds.convert(f'{shortname}.xyz', fmt='extxyz')
  """
}
