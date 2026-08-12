import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// matplotlib viridis, linearly interpolated; t in [0,1] -> [r,g,b] in [0,1]
const VIRIDIS = [
  [0.267, 0.005, 0.329], [0.283, 0.141, 0.458], [0.254, 0.265, 0.530], [0.207, 0.372, 0.553],
  [0.164, 0.471, 0.558], [0.128, 0.567, 0.551], [0.135, 0.659, 0.518], [0.267, 0.749, 0.441],
  [0.478, 0.821, 0.318], [0.741, 0.873, 0.150], [0.993, 0.906, 0.144],
];
export function viridis(t) {
  t = Math.max(0, Math.min(1, t));
  const i = Math.min(9, Math.floor(t * 10)), f = t * 10 - i;
  const a = VIRIDIS[i], b = VIRIDIS[i + 1];
  return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
}

// Canvas sized to the text so long labels are never cut off.
function textCanvas(text) {
  const c = document.createElement('canvas');
  const g = c.getContext('2d');
  const font = '28px ui-monospace, monospace';
  g.font = font;
  c.width = Math.ceil(g.measureText(text).width) + 16;
  c.height = 64;
  g.font = font; // resizing the canvas resets context state
  g.fillStyle = '#9aa7b5';
  g.textAlign = 'center';
  g.textBaseline = 'middle';
  g.fillText(text, c.width / 2, c.height / 2);
  return c;
}

// Billboard label (always faces the camera). `size` = world-space text height.
export function textSprite(text, size = 0.5) {
  const c = textCanvas(text);
  const s = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(c), transparent: true }));
  s.scale.set(size * c.width / c.height, size, 1);
  return s;
}
// Shared renderer/camera/controls shell for all views.
// - orthographic camera by default, frustum matched to the perspective view (locked decision, map #1)
// - continuous render loop with optional per-frame tick (views use it for play/marker animation)
export class SceneBase {
  constructor(container, { camPos = [7, 5, 8], target = [0, 0, 0], ortho = true } = {}) {
    this.container = container;
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x0b0e13);

    const aspect = container.clientWidth / Math.max(1, container.clientHeight);
    this.persp = new THREE.PerspectiveCamera(50, aspect, 0.1, 300);
    this.persp.position.set(...camPos);
    this.targetV = new THREE.Vector3(...target);
    this.orthoHalfH = this.persp.position.distanceTo(this.targetV) * Math.tan(this.persp.fov * Math.PI / 360);
    this.orthoCam = new THREE.OrthographicCamera(
      -this.orthoHalfH * aspect, this.orthoHalfH * aspect, this.orthoHalfH, -this.orthoHalfH, -300, 300);
    this.orthoCam.position.set(...camPos);
    this.camera = ortho ? this.orthoCam : this.persp;

    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setPixelRatio(window.devicePixelRatio);
    this.renderer.setSize(container.clientWidth, container.clientHeight);
    container.appendChild(this.renderer.domElement);

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.target.copy(this.targetV);

    this.scene.add(new THREE.AmbientLight(0xffffff, 0.7));
    const dir = new THREE.DirectionalLight(0xffffff, 1.2);
    dir.position.set(5, 10, 3);
    this.scene.add(dir);

    this.onTick = null;
    this._clock = new THREE.Clock();
    this.renderer.setAnimationLoop(() => {
      const dt = this._clock.getDelta();
      if (this.onTick) this.onTick(dt);
      this.controls.update();
      this.renderer.render(this.scene, this.camera);
    });

    this._ro = new ResizeObserver(() => this._resize());
    this._ro.observe(container);
    this._ray = new THREE.Raycaster();
    this._mouse = new THREE.Vector2();
  }

  _resize() {
    const w = this.container.clientWidth, h = Math.max(1, this.container.clientHeight);
    const aspect = w / h;
    this.persp.aspect = aspect;
    this.persp.updateProjectionMatrix();
    this.orthoCam.left = -this.orthoHalfH * aspect;
    this.orthoCam.right = this.orthoHalfH * aspect;
    this.orthoCam.updateProjectionMatrix();
    this.renderer.setSize(w, h);
  }

  setOrtho(on) {
    const next = on ? this.orthoCam : this.persp;
    if (next === this.camera) return;
    next.position.copy(this.camera.position);
    this.camera = next;
    this.camera.updateProjectionMatrix();
    this.controls.object = this.camera;
    this.controls.update();
  }

  // Snap the camera to an axis-aligned projection: 'yx' | 'zx' | 'zy'
  // (first axis = screen up, second = screen right), keeping the current distance and target.
  snapTo(plane) {
    const P = {
      yx: { pos: [0, 0, 1], up: [0, 1, 0] },
      zx: { pos: [0, -1, 0], up: [0, 0, 1] },
      zy: { pos: [1, 0, 0], up: [0, 0, 1] },
    }[plane];
    if (!P) return;
    const t = this.controls.target;
    const dist = this.camera.position.distanceTo(t);
    const pos = new THREE.Vector3(...P.pos).multiplyScalar(dist).add(t);
    for (const cam of [this.persp, this.orthoCam]) {
      cam.up.set(...P.up);
      cam.position.copy(pos);
      cam.lookAt(t);
    }
    this.controls.update();
  }

  // Raycast a pointer event against an object; returns the intersection or null.
  raycast(event, object) {
    const r = this.renderer.domElement.getBoundingClientRect();
    this._mouse.set(((event.clientX - r.left) / r.width) * 2 - 1, -((event.clientY - r.top) / r.height) * 2 + 1);
    this._ray.setFromCamera(this._mouse, this.camera);
    return this._ray.intersectObject(object)[0] ?? null;
  }

  dispose() {
    this._ro.disconnect();
    this.renderer.setAnimationLoop(null);
    this.renderer.dispose();
    this.scene.traverse((o) => {
      o.geometry?.dispose?.();
      if (o.material) (Array.isArray(o.material) ? o.material : [o.material]).forEach((m) => { m.map?.dispose?.(); m.dispose?.(); });
    });
    this.renderer.domElement.remove();
  }
}
