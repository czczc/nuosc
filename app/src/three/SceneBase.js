import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { theme } from '../theme.js';

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

// Global color palettes (store.palette), t in [0,1] -> [r,g,b] in [0,1].
const clamp01 = (v) => Math.max(0, Math.min(1, v));

// classic jet-style rainbow: blue (low) -> cyan -> green -> yellow -> red (high)
function rainbow(t) {
  t = clamp01(t);
  return [clamp01(1.5 - Math.abs(4 * t - 3)), clamp01(1.5 - Math.abs(4 * t - 2)), clamp01(1.5 - Math.abs(4 * t - 1))];
}
// diverging cool-warm: blue (low) -> light gray -> red (high)
function coolwarm(t) {
  t = clamp01(t);
  const a = [0.23, 0.30, 0.75], m = [0.86, 0.86, 0.86], b = [0.71, 0.02, 0.15];
  const [p, q, f] = t < 0.5 ? [a, m, 2 * t] : [m, b, 2 * t - 1];
  return [p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1]), p[2] + f * (q[2] - p[2])];
}
function grayscale(t) {
  t = clamp01(t);
  return [t, t, t];
}
export const PALETTES = { rainbow, viridis, coolwarm, grayscale };

// Canvas sized to the text so long labels are never cut off.
// text is a string or an array of [text, color?] segments (multi-color labels).
function textCanvas(text, color) {
  const segs = Array.isArray(text) ? text : [[text, color]];
  const c = document.createElement('canvas');
  const g = c.getContext('2d');
  const font = '28px ui-monospace, monospace';
  g.font = font;
  const widths = segs.map(([t]) => g.measureText(t).width);
  c.width = Math.ceil(widths.reduce((a, b) => a + b, 0)) + 16;
  c.height = 64;
  g.font = font; // resizing the canvas resets context state
  g.textAlign = 'left';
  g.textBaseline = 'middle';
  let x = 8;
  for (let i = 0; i < segs.length; i++) {
    g.fillStyle = segs[i][1] ?? theme().label;
    g.fillText(segs[i][0], x, c.height / 2);
    x += widths[i];
  }
  return c;
}

// Billboard label (always faces the camera). SceneBase rescales it every frame so
// its on-screen height stays constant while zooming; `k` scales relative to the
// standard label height.
export function textSprite(text, k = 1, color) {
  const c = textCanvas(text, color);
  const s = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(c), transparent: true }));
  s.userData.labelK = k;
  s.userData.labelAspect = c.width / c.height;
  return s;
}

// Label height as a fraction of the viewport height.
const LABEL_FRAC = 0.034;
// Shared renderer/camera/controls shell for all views.
// - always orthographic, frustum matched to a reference perspective view (locked decision, map #1)
// - continuous render loop with optional per-frame tick (views use it for play/marker animation)
export class SceneBase {
  constructor(container, { camPos = [7, 5, 8], target = [0, 0, 0], snapPlanes = {} } = {}) {
    this.container = container;
    this._snapPlanes = { front: [0, 0, 1], top: [0, 1, 0.001], side: [-1, 0, 0], ...snapPlanes };
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(theme().stage);

    const aspect = container.clientWidth / Math.max(1, container.clientHeight);
    this.persp = new THREE.PerspectiveCamera(50, aspect, 0.1, 300);
    this.persp.position.set(...camPos);
    this.targetV = new THREE.Vector3(...target);
    this.orthoHalfH = this.persp.position.distanceTo(this.targetV) * Math.tan(this.persp.fov * Math.PI / 360);
    this.orthoCam = new THREE.OrthographicCamera(
      -this.orthoHalfH * aspect, this.orthoHalfH * aspect, this.orthoHalfH, -this.orthoHalfH, -300, 300);
    this.orthoCam.position.set(...camPos);
    this.camera = this.orthoCam;

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
    this._snapAnim = null;
    this.controls.addEventListener('start', () => { this._snapAnim = null; });
    this._clock = new THREE.Clock();
    this.renderer.setAnimationLoop(() => {
      const dt = this._clock.getDelta();
      if (this.onTick) this.onTick(dt);
      if (this._snapAnim) this._stepSnap(dt);
      else this.controls.update();
      this._updateLabels();
      this.renderer.render(this.scene, this.camera);
    });

    this._ro = new ResizeObserver(() => this._resize());
    this._ro.observe(container);
    this._ray = new THREE.Raycaster();
    this._mouse = new THREE.Vector2();
    this._tmpV = new THREE.Vector3();
  }

  // Rescale label sprites so their on-screen height stays LABEL_FRAC of the
  // viewport regardless of camera zoom/distance.
  _updateLabels() {
    const ortho = this.camera === this.orthoCam;
    const orthoH = (this.orthoCam.top - this.orthoCam.bottom) / this.orthoCam.zoom;
    const perspK = 2 * Math.tan(this.persp.fov * Math.PI / 360);
    this.scene.traverse((o) => {
      if (!o.isSprite || !o.userData.labelK) return;
      const visH = ortho ? orthoH : perspK * this.persp.position.distanceTo(o.getWorldPosition(this._tmpV));
      const s = LABEL_FRAC * visH * o.userData.labelK;
      o.scale.set(s * o.userData.labelAspect, s, 1);
    });
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

  // Snap the camera to an axis-aligned projection, keeping the current distance and target.
  // Defaults: 'front' from +z (x right, y up), 'top' from +y (x right, -z up),
  // 'side' from -x (z right, y up); views override per-plane via the snapPlanes option.
  // Camera up stays world +y — OrbitControls caches its up axis at construction, so
  // changing it makes orbiting feel twisted. The top view is tilted by an epsilon
  // instead, so lookAt/OrbitControls never degenerate.
  snapTo(plane) {
    const to = this._snapPlanes[plane];
    if (!to) return;
    const t = this.controls.target;
    const dist = this.camera.position.distanceTo(t);
    const endPos = new THREE.Vector3(...to).normalize().multiplyScalar(dist).add(t);
    const m = new THREE.Matrix4().lookAt(endPos, t, new THREE.Vector3(0, 1, 0));
    this._snapAnim = {
      q0: this.camera.quaternion.clone(),
      q1: new THREE.Quaternion().setFromRotationMatrix(m),
      dist,
      u: 0,
    };
  }

  // Advance the snap tween: slerp the full camera orientation (roll included —
  // a position-only tween jumps at the top view, where lookAt's roll is degenerate)
  // and place the camera on its -z axis at constant distance from the target,
  // easing with smoothstep. Cancelled when the user grabs the controls; the loop
  // skips controls.update() while active so it can't overwrite the pose.
  _stepSnap(dt) {
    const a = this._snapAnim;
    if (!a) return;
    a.u = Math.min(1, a.u + dt / 0.45);
    const k = a.u * a.u * (3 - 2 * a.u);
    const q = new THREE.Quaternion().slerpQuaternions(a.q0, a.q1, k);
    const pos = new THREE.Vector3(0, 0, 1).applyQuaternion(q).multiplyScalar(a.dist).add(this.controls.target);
    for (const cam of [this.persp, this.orthoCam]) {
      cam.up.set(0, 1, 0);
      cam.position.copy(pos);
      cam.quaternion.copy(q);
    }
    if (a.u >= 1) this._snapAnim = null;
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
