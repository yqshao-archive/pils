// this workflow generates the input geometries for systems with different sizes
// Naming of the system: ${acid}-${base}-a${nacid}b${nbase}i${nion}-rho${density}
nextflow.enable.dsl=2

params.publish = 'build'
params.acid = 'hoac'
params.base = 'c1im'

process moltemplate {
  tag "$name"
  publishDir "$params.publish/$name"
  label 'moltemplate'

  input:
    val name
    tuple val(acid), val(base), val(anion), val(cation), \
          val(nacid), val(nbase), val(nion), val(box), path(files)

  output:
    tuple val(name), path("system.{in.init,in.settings,data}")

  """
  packmol << EOF > packmol.log
    tolerance  2.0
    filetype   pdb
    output     system.pdb

  ${nacid!=0 ? """
  structure  ${acid}.pdb
    number   ${nacid}
    inside   cube  0. 0. 0. ${box-2.0}
  end structure""" : ""}

  ${nbase!=0 ? """
  structure  ${base}.pdb
    number   ${nbase}
    inside   cube  0. 0. 0. ${box-2.0}
  end structure""" : ""}

  ${nion!=0 ? """
  structure  ${anion}.pdb
    number   ${nion}
    inside   cube  0. 0. 0. ${box-2.0}
  end structure""" : ""}

  ${nion!=0 ? """
  structure  ${cation}.pdb
    number   ${nion}
    inside   cube  0. 0. 0. ${box-2.0}
  end structure""" : ""}

  EOF

  #
  cat << EOF > system.lt
  import "${acid}.lt"
  import "${base}.lt"
  import "${anion}.lt"
  import "${cation}.lt"

  # Periodic boundary conditions:
  write_once("Data Boundary") {
     0.0  $box  xlo xhi
     0.0  $box  ylo yhi
     0.0  $box  zlo zhi
  }

  # Acids (if any)
  ${nacid!=0 ? "acids = new $acid[$nacid]" : ""}

  # Bases (if any)
  ${nbase!=0 ? "bases = new $base[$nbase]" : ""}

  # Ions (if any)
  ${nion!=0 ? "anions = new $anion[$nion]" : ""}
  ${nion!=0 ? "cations = new $cation[$nion]" : ""}

  EOF

  moltemplate.sh -pdb system.pdb system.lt
  cleanup_moltemplate.sh
  """
}

def tag2inp (tag, rho) {
  // converting the strings to variables
  def acid = params.acid
  def base = params.base
  def anion = acid[1..-1]
  def cation = 'h'+base
  def matches = (tag =~ /a(\d+)b(\d+)i(\d+)/ )[0]
  def nacid = matches[1] as Integer
  def nbase = matches[2] as Integer
  def nion  = matches[3] as Integer
  def name = "$acid-$base-$tag-rho$rho"
  rho = rho as BigDecimal // g cm^-3

  // compute density; molecular weights are hard coded here
  def masses = ['hoac': 60.052, 'c1im': 82.108] // g mol^-1
  def NA = 6.022e23 // mol^-1
  def mass = masses[acid]*(nion+nacid) + masses[base]*(nion+nbase)
  def box = Math.cbrt(mass*1e24/rho/NA).round(6) // cm -> Ang

  // input files
  def files = file("skel/template/{${[acid,base,cation,anion].join(',')}}.{lt,pdb}")
  return [name, [acid, base, anion, cation, nacid, nbase, nion, box, files]]
}
