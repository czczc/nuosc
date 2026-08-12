// Minimal matplotlib-style scaffold for the 2D companion canvases.
// plot2d() draws the frame, grid, ticks and outside axis titles, and returns
// data->pixel mappers plus the plot rectangle for the caller's curves.

function niceStep(span, n) {
  const raw = span / n;
  const mag = 10 ** Math.floor(Math.log10(raw));
  for (const m of [1, 2, 2.5, 5, 10]) if (m * mag >= raw) return m * mag;
  return 10 * mag;
}

function ticks(a, b, n) {
  const step = niceStep(b - a, n);
  const out = [];
  for (let t = Math.ceil(a / step - 1e-9) * step; t <= b + 1e-9 * step; t += step) out.push(t);
  return { at: out, step };
}

function fmtTick(t, step) {
  if (t === 0) return '0';
  let prec = 0;
  while (prec < 6 && Math.abs(Math.round(step * 10 ** prec) - step * 10 ** prec) > 1e-6) prec++;
  return t.toFixed(prec);
}

export function plot2d(ctx, w, h, dpr, { x: [x0, x1], y: [y0, y1], xTitle, yTitle, xTicks = 6, yTicks = 4 }) {
  const ml = 48 * dpr, mr = 10 * dpr, mt = 8 * dpr, mb = 30 * dpr;
  const pw = w - ml - mr, ph = h - mt - mb;
  const X = (x) => ml + pw * (x - x0) / (x1 - x0);
  const Y = (y) => mt + ph * (1 - (y - y0) / (y1 - y0));

  const gridCol = 'rgba(140,155,170,0.15)', frameCol = 'rgba(140,155,170,0.5)', ink = '#9aa7b5';
  const xt = ticks(x0, x1, xTicks), yt = ticks(y0, y1, yTicks);
  ctx.lineWidth = dpr;

  ctx.strokeStyle = gridCol;
  for (const t of xt.at) { const p = X(t); ctx.beginPath(); ctx.moveTo(p, mt); ctx.lineTo(p, mt + ph); ctx.stroke(); }
  for (const t of yt.at) { const p = Y(t); ctx.beginPath(); ctx.moveTo(ml, p); ctx.lineTo(ml + pw, p); ctx.stroke(); }

  const frame = () => { ctx.strokeStyle = frameCol; ctx.lineWidth = dpr; ctx.strokeRect(ml, mt, pw, ph); };
  frame();

  ctx.font = `${9 * dpr}px ui-monospace, monospace`;
  ctx.fillStyle = ink;
  ctx.strokeStyle = frameCol;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  for (const t of xt.at) {
    const p = X(t);
    ctx.beginPath(); ctx.moveTo(p, mt + ph); ctx.lineTo(p, mt + ph + 3 * dpr); ctx.stroke();
    ctx.fillText(fmtTick(t, xt.step), p, mt + ph + 5 * dpr);
  }
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  for (const t of yt.at) {
    const p = Y(t);
    ctx.beginPath(); ctx.moveTo(ml - 3 * dpr, p); ctx.lineTo(ml, p); ctx.stroke();
    ctx.fillText(fmtTick(t, yt.step), ml - 5 * dpr, p);
  }

  ctx.font = `${10 * dpr}px ui-monospace, monospace`;
  ctx.textAlign = 'center';
  if (xTitle) { ctx.textBaseline = 'alphabetic'; ctx.fillText(xTitle, ml + pw / 2, h - 5 * dpr); }
  if (yTitle) {
    ctx.save();
    ctx.translate(11 * dpr, mt + ph / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = 'middle';
    ctx.fillText(yTitle, 0, 0);
    ctx.restore();
  }
  ctx.textAlign = 'left';
  ctx.textBaseline = 'alphabetic';

  return { X, Y, ml, mt, pw, ph, frame };
}

// Colored legend, top-right inside the plot, on a dark backdrop chip.
export function legend(ctx, dpr, P, items) {
  ctx.font = `${10 * dpr}px ui-monospace, monospace`;
  const wds = items.map((it) => ctx.measureText(it.label).width);
  const total = wds.reduce((a, b) => a + b, 0) + 10 * dpr * (items.length - 1);
  let x = P.ml + P.pw - total - 8 * dpr;
  ctx.fillStyle = 'rgba(11,14,19,0.75)';
  ctx.fillRect(x - 4 * dpr, P.mt + 3 * dpr, total + 8 * dpr, 14 * dpr);
  for (let i = 0; i < items.length; i++) {
    ctx.fillStyle = items[i].color;
    ctx.fillText(items[i].label, x, P.mt + 14 * dpr);
    x += wds[i] + 10 * dpr;
  }
}
