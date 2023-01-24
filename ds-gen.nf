#!/usr/bin/env nextflow
nextflow.enable.dsl=2
nextflow.preview.recursion=true
/*
 This is the workflow for an trial run, it has the following entry points
 `-entry lmpinit` generates the initial geometry guess with moltemplate and lammps
 `-entry cp2k` generates the cp2k trajectories or extends the trajectory
 `-entry train` trains the inital MD models and perform some basic analyses
*/

//======================================//
// THE GENERATION OF INITIAL GEOMETRIES //
//======================================//
include { moltemplate; tag2inp } from './nf/build.nf' addParams(publish: 'trajs/build')
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
params.cp2k_tag = 'nvt-scan'
params.cp2k_aux = 'skel/cp2k-aux/*'
params.cp2k_geo = 'trajs/lmp/*/equi.dump'
params.cp2k_each_ps  = '2'
params.cp2k_max_ps = '100'
include { cp2kGenInp; cp2k } from './tips/nextflow/cp2k.nf' addParams(publish: 'trajs/cp2k')


def name2ps (name) {(name=~/-(\d+)ps/)[0][1].toInteger()}
def name2geo (name) {(name=~/-\d+ps\/(.*)/)[0][1]}

workflow cp2kprod {
  each_ps = params.cp2k_each_ps.toInteger()
  max_ps = params.cp2k_max_ps.toInteger()
  channel.fromPath(params.cp2k_geo) \
    | map {it -> ["$params.cp2k_tag-0-${each_ps}ps/${it.parent.name}",
                  file("skel/cp2k/$params.cp2k_tag-${each_ps}ps.inp"),
                  it,
                  "--fmt lammps-dump --idx -1 --emap ${file("trajs/build/${it.parent.name}/system.data")}"]} \
    | cp2kGenInp \
    | collect(flat:false) \
    | set {ch_cp2k_init}

  cp2kloop.recurse(ch_cp2k_init).until{name2ps(it)>max_ps}
}

workflow cp2kloop {
  take: ch_inp

  main:
  each_ps = params.cp2k_each_ps.toInteger()
  cp2k_aux = file(params.cp2k_aux)
  inp_size = file(params.cp2k_geo).size
  ch_inp \
    | flatMap {inps -> \
               inps.collect{ name, inp -> \
                            [name, inp, cp2k_aux]}} \
    | cp2k

  cp2k.out.restart \
    | map {name, files -> \
           ["$params.cp2k_tag-${name2ps(name)}-${name2ps(name)+each_ps}-ps/${name2geo(name)}",
            files.sort {f -> (f=~/_(\d+).restart/)[0][1].toInteger()}]} \
    | collate(inp_size) \
    | set {nx_inp}

  emit:
  nx_inp
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
