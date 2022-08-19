// PiView.js
//
// This is a visualization module to generate simple visualizations of p(roperties)
// and i(nteractions). This is created for the purpose of visualizting the
// AML predictions in the manuscript. A updated version might be published later
// separately.
//
// written by Yunqi Shao 2021-11-29
import * as THREE from "https://cdn.skypack.dev/three@0.135.0";
import { TrackballControls } from "https://cdn.skypack.dev/three@0.135.0/examples/jsm/controls/TrackballControls.js"
import { RoomEnvironment } from "https://cdn.skypack.dev/three@0.135.0/examples/jsm/environments/RoomEnvironment.js"

var radii = {'vdW':
             [1.2, 1.2, 1.4, 1.82, 1.53, 1.92, 1.7, 1.55, 1.52, 1.47, 1.54, 2.27, 1.73, 1.84, 2.1, 1.8, 1.8, 1.75, 1.88],
             'covalent':
             [0.2, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06]};

var jmol = ["#ff0000", "#ffffff", "#d9ffff", "#cc80ff", "#c2ff00", "#ffb5b5", "#909090", "#3050f8", "#ff0d0d", "#90e050", "#b3e3f5", "#ab5cf2", "#8aff00", "#bfa6a6", "#f0c8a0", "#ff8000", "#ffff30", "#1ff01f", "#80d1e3", "#8f40d4", "#3dff00", "#e6e6e6", "#bfc2c7"];
var trans     = new THREE.MeshBasicMaterial({color: 0xffffff, side: THREE.DoubleSide, transparent: true, opacity: 0 });
var black     = new THREE.MeshBasicMaterial({color: 0x000000, side: THREE.BackSide});
var materials = {};
var zoom      = zoom      ?? 1;
var margin    = margin    ?? 0.1;
var quality   = quality   ?? 24;  // subdivision for balls and bonds
var atomScale = atomScale ?? 0.24; //
var edge      = edge      ?? 0.1;
var bondRadi  = bondRadi  ?? 0.15;
var radiusType= radiusType?? 'vdW';
var luxurious = luxurious ?? false;
var radi = radii[radiusType];

if (luxurious){edge=0.0};

function get_mat(color){
  if (!materials[color]) {
      materials[color] = new THREE.MeshBasicMaterial({color: color});
      if (luxurious) {
          materials[color] = new THREE.MeshStandardMaterial({color: color, roughness: 0,metalness: 1})}
  };
  return materials[color]
}

function animate(renderer, scene, camera, controls) {
    requestAnimationFrame(()=> animate(renderer, scene, camera, controls) );
    controls.update();
    renderer.render(scene, camera);
}

function renderPiView (node) {
    // canvas setup
    const canvas  = document.createElement('canvas');
    canvas.style  = node.getAttribute('cstyle')?? "border: 0.2em solid";
    canvas.width  = node.getAttribute('width') ?? 200;
    canvas.height = node.getAttribute('height')?? 200;
    // load atoms first compute the default camera limits
    const atomsDict = JSON.parse(node.getAttribute('atoms'));
    var wlim = Math.max(...atomsDict.elems.map((e,i)=>{
        return Math.abs(atomsDict.coord[i][0])+radi[e]*atomScale}))+edge+margin;
    var hlim = Math.max(...atomsDict.elems.map((e,i)=>{
        return Math.abs(atomsDict.coord[i][1])+radi[e]*atomScale}))+edge+margin;
    wlim = Math.max(wlim, hlim/canvas.height*canvas.width)/zoom;
    hlim = Math.max(hlim, wlim/canvas.width*canvas.height)/zoom;
    // three.js setup
    const camera   = new THREE.OrthographicCamera(-wlim, wlim, hlim, -hlim, 1, 2000);
    const renderer = new THREE.WebGLRenderer({canvas: canvas, antialias: true, alpha: true,});
    const scene    = new THREE.Scene();
    camera.position.z = 1000;
    // make atoms
    var tmp = makeAtoms(atomsDict);
    var atoms = tmp[0];
    var outline = tmp[1];

    if (luxurious) {
        var pmremGenerator = new THREE.PMREMGenerator( renderer );
        var environment = pmremGenerator.fromScene( new RoomEnvironment(), 0.04 ).texture;
        scene.environment = environment;}

    const props = JSON.parse(node.getAttribute('props'));
    if (props){
        atoms.children.forEach(
            (atom, idx)=>{atom.material=get_mat(props[idx])})};

    scene.add(atoms);
    scene.add(outline);
    const inter = JSON.parse(node.getAttribute('inter'));
    if (inter){
        var bond_group = new THREE.Group();
        inter.forEach(
            (bond, idx)=>{
                makeBond(atomsDict.coord[bond[0]],
                         atomsDict.coord[bond[1]],
                         Math.sqrt((radi[atomsDict.elems[bond[0]]]*atomScale+edge)**2-bondRadi**2),
                         Math.sqrt((radi[atomsDict.elems[bond[1]]]*atomScale+edge)**2-bondRadi**2),
                         bond[2], bond_group)})
        scene.add(bond_group)};


    // // render
    // renderer.setClearColor(0xffffff);
    var download = node.getAttribute('download') ?? false;
    if (download) {
        screenshot(renderer, scene, camera, download);
    }else{
        node.appendChild(canvas);
        const controls = new TrackballControls(camera, canvas);
        controls.noPan = true;
        controls.zoomSpeed = 1;
        const button  = document.createElement('button');
        button.style  = "z-index:2; position: absolute; top:0px; left:0px;";
        button.innerHTML  = 'Shot';
        button.addEventListener('click', ()=>screenshot(renderer, scene, camera, null));
        //node.appendChild(button);
        animate(renderer, scene, camera, controls)}}

function makeAtom(elem, coord, group, out_group) {
    const geo  = new THREE.SphereGeometry(radi[elem]*atomScale, quality, quality);
    const mat  = get_mat(jmol[elem]);
    const ball = new THREE.Mesh(geo, mat);
    const out_geo  = new THREE.SphereGeometry(radi[elem]*atomScale + edge, quality, quality);
    const out_mat  = get_mat(0x000000);
    const out_ball = new THREE.Mesh(out_geo, black);
    ball.position.set(...coord);
    out_ball.position.set(...coord);
    group.add(ball);
    out_group.add(out_ball)}

function makeBond(coord1, coord2, offset1, offset2, color, group) {
    var coord1 = new THREE.Vector3(...coord1);
    var coord2 = new THREE.Vector3(...coord2);
    var diff = new THREE.Vector3();
    diff.subVectors(coord2, coord1);
    var distance = diff.length();
    coord1.addScaledVector(diff, offset1/distance)
    diff.addScaledVector(diff, -offset2/distance-offset1/distance)
    diff.addScaledVector(diff, -0.5)
    distance=diff.length();
    var HALF_PI = Math.PI * 0.5;
    var position = coord1.add(diff.divideScalar(2));
    var geo = new THREE.CylinderGeometry(bondRadi, bondRadi, distance, quality, 1, true);
    var orientation = new THREE.Matrix4();
    var offsetRotation = new THREE.Matrix4();
    var offsetPosition = new THREE.Matrix4();
    orientation.lookAt(new THREE.Vector3(0, 0, 0), diff, new THREE.Vector3(0, 1, 0));
    offsetRotation.makeRotationX(HALF_PI);
    orientation.multiply(offsetRotation);
    geo.applyMatrix4(orientation);
    var bond = new THREE.Mesh(geo, [get_mat(color), trans, get_mat(color)]);
    bond.position.set(position.x, position.y, position.z);
    group.add(bond)}

function makeAtoms(atoms){
  var atom_group = new THREE.Group();
  var out_group = new THREE.Group();
  atoms.elems.map((elem, i) => makeAtom(elem, atoms.coord[i], atom_group, out_group))
  atom_group.matrixAutoUpdate = false;
  out_group.matrixAutoUpdate  = false;
  return [atom_group, out_group]}

function screenshot(renderer, scene, camera, name) {
    var name = name ?? 'screenshot'
    var a = document.createElement('a');
    renderer.render(scene, camera);
    a.href = renderer.domElement.toDataURL().replace("image/png", "image/octet-stream");
    a.download = name+'.png'
    a.click()}

var node = document.getElementById(uuid);
renderPiView(node)
