// --- Physics Constants (Fundamental) ---
const DEG_TO_RAD = Math.PI / 180.0;
const KM_TO_GEV_INV = 5.0677e9; // Conversion: 1 km = 5.0677e9 GeV^-1 (using hbar*c = 1)
const GF = 1.1663787e-5; // Fermi constant in GeV^-2
const A_FACTOR_EV2 = 0.76e-4; // eV² / (g/cm³ * GeV)
const YE = 0.5; // Electron fraction
const CONV_FACTOR_PHASE = 1.267 * 2; // Note: Factor of 2 seems unusual here, typical is 1.267. Keeping as is for now.

// --- DOM Elements ---
// Existing Controls
const dcpSlider = document.getElementById('dcpSlider');
const dcpInput = document.getElementById('dcpInput');
const dcpValueSpan = document.getElementById('dcpValue');
// Mass Ordering Toggle
const massOrderToggle = document.getElementById('massOrderToggle');
const massOrderLabel = document.getElementById('massOrderLabel');
// Other Controls
const distanceSlider = document.getElementById('distanceSlider');
const distanceInput = document.getElementById('distanceInput');
const distanceValueSpan = document.getElementById('distanceValue');
const densitySlider = document.getElementById('densitySlider');
const densityInput = document.getElementById('densityInput');
const densityValueSpan = document.getElementById('densityValue');
// New Physics Parameter Controls
const theta12Slider = document.getElementById('theta12Slider');
const theta12Input = document.getElementById('theta12Input');
const theta12ValueSpan = document.getElementById('theta12Value');
const theta13Slider = document.getElementById('theta13Slider');
const theta13Input = document.getElementById('theta13Input');
const theta13ValueSpan = document.getElementById('theta13Value');
const theta23Slider = document.getElementById('theta23Slider');
const theta23Input = document.getElementById('theta23Input');
const theta23ValueSpan = document.getElementById('theta23Value');
const dm2_21Slider = document.getElementById('dm2_21Slider');
const dm2_21Input = document.getElementById('dm2_21Input');
const dm2_21ValueSpan = document.getElementById('dm2_21Value');
const dm2_31Slider = document.getElementById('dm2_31Slider');
const dm2_31Input = document.getElementById('dm2_31Input');
const dm2_31ValueSpan = document.getElementById('dm2_31Value');
// Chart Context
const ctx = document.getElementById('oscillationChart').getContext('2d');

// --- Global Variables for Physics Parameters & Derived Values ---
let currentTheta12Rad, currentTheta13Rad, currentTheta23Rad;
let currentS12_2, currentC12_2, currentS13_2, currentC13_2, currentS23_2, currentC23_2;
let currentSin2_2T13, currentSin2_2T12, currentSin2_2T23;
let currentDm2_21, currentDm2_31, currentDm2_32;
let currentSin2T12, currentSin2T13, currentSin2T23; // Need sin(2*theta) terms as well

// --- Function to Update Physics Parameters from Inputs ---
function updatePhysicsParameters() {
    // Read degrees from inputs
    const theta12Deg = parseFloat(theta12Input.value);
    const theta13Deg = parseFloat(theta13Input.value);
    const theta23Deg = parseFloat(theta23Input.value);

    // Convert angles to radians
    currentTheta12Rad = theta12Deg * DEG_TO_RAD;
    currentTheta13Rad = theta13Deg * DEG_TO_RAD;
    currentTheta23Rad = theta23Deg * DEG_TO_RAD;

    // Pre-calculate sin/cos terms
    const s12 = Math.sin(currentTheta12Rad); const c12 = Math.cos(currentTheta12Rad);
    const s13 = Math.sin(currentTheta13Rad); const c13 = Math.cos(currentTheta13Rad);
    const s23 = Math.sin(currentTheta23Rad); const c23 = Math.cos(currentTheta23Rad);

    // Pre-calculate squared terms
    currentS12_2 = s12**2; currentC12_2 = c12**2;
    currentS13_2 = s13**2; currentC13_2 = c13**2;
    currentS23_2 = s23**2; currentC23_2 = c23**2;

    // Pre-calculate combinations for formula
    currentSin2T12 = Math.sin(2 * currentTheta12Rad); // Needed for interference term
    currentSin2T13 = Math.sin(2 * currentTheta13Rad); // Needed for interference term
    currentSin2T23 = Math.sin(2 * currentTheta23Rad); // Needed for interference term
    currentSin2_2T13 = currentSin2T13**2;
    currentSin2_2T12 = currentSin2T12**2;
    currentSin2_2T23 = currentSin2T23**2;

    // Read mass splittings magnitude (and unscale)
    const dm2_21_mag = parseFloat(dm2_21Input.value) * 1e-5;
    const dm2_31_mag = parseFloat(dm2_31Input.value) * 1e-3;

    // Set sign based on mass ordering toggle
    const isNormalOrdering = massOrderToggle.checked;
    currentDm2_21 = dm2_21_mag; // Always positive
    currentDm2_31 = isNormalOrdering ? dm2_31_mag : -dm2_31_mag;

    // Calculate dm2_32 based on the signed dm2_31
    currentDm2_32 = currentDm2_31 - currentDm2_21;

    // Update toggle label and color
    if (isNormalOrdering) {
        massOrderLabel.textContent = 'Normal';
        massOrderLabel.style.color = ''; // Reset to default CSS color
    } else {
        massOrderLabel.textContent = 'Inverted';
        massOrderLabel.style.color = '#656060'; // Set specific color for Inverted
    }
}

// --- Chart Setup ---
let oscillationChart; // To hold the Chart object

// Helper function for sinc(x) = sin(x)/x
function sinc(x) {
    if (Math.abs(x) < 1e-9) {
        return 1.0 - x*x/6.0 + x*x*x*x/120.0; // Taylor expansion for small x
    } else {
        return Math.sin(x) / x;
    }
}

// ** Function to calculate oscillation probability using dynamically updated parameters **
function calculateOscillationProbabilityMatter(energyGeV, distanceKm, density_g_cm3, deltaCpDeg, isAntineutrino = false) {
    if (energyGeV <= 0) return 0;

    // 1. Calculate Matter Potential A in eV² (flip sign for antineutrinos)
    const matterSign = isAntineutrino ? -1.0 : 1.0;
    const A_eV2 = matterSign * A_FACTOR_EV2 * density_g_cm3 * YE * energyGeV;

    // 2. Calculate vacuum oscillation phases / arguments using currentDm2_31
    const phase_31 = CONV_FACTOR_PHASE * currentDm2_31 * distanceKm / energyGeV;
    const D31 = phase_31 / 2.0; // Δ = Δm²₍₃₁₎ L / (4E)

    // 3. Define Ahat and α using current mass splittings
    const Ahat = (Math.abs(currentDm2_31) < 1e-20) ? 0 : A_eV2 / currentDm2_31;
    const alpha = (Math.abs(currentDm2_31) < 1e-20) ? 0 : currentDm2_21 / currentDm2_31;
    const deltaCpRad = deltaCpDeg * DEG_TO_RAD;
    // For neutrinos: interference phase = D31 + δCP, for antineutrinos: D31 − δCP
    const interferencePhase = D31 + (isAntineutrino ? -deltaCpRad : deltaCpRad);

    // Calculate arguments for the sine factors
    const arg_1_minus_Ahat = D31 * (1.0 - Ahat);
    const arg_Ahat = D31 * Ahat;

    // Factor for the atmospheric term (sin((1-Ahat)Δ)/(1-Ahat))
    const factor1 = (Math.abs(1.0 - Ahat) < 1e-9) ? D31 : Math.sin(arg_1_minus_Ahat) / (1.0 - Ahat);
    // Use current S23_2 and currentSin2_2T13
    const term1_pdg = currentS23_2 * currentSin2_2T13 * factor1 * factor1;

    // Factor for the solar term (sin(AhatΔ)/Ahat)
    const factor2 = (Math.abs(Ahat) < 1e-9) ? D31 : Math.sin(arg_Ahat) / Ahat;
    // Use current C23_2 (derived from S23_2) and currentSin2_2T12
    const term2_pdg = alpha * alpha * (1.0 - currentS23_2) * currentSin2_2T12 * factor2 * factor2; // C23_2 = 1 - S23_2

    // Interference term using the corrected phase: cos(D31 ± δCP)
    // Use current sin(2*theta) terms
    const term3_pdg = alpha * currentSin2T12 * currentSin2T13 * currentSin2T23 *
                        Math.cos(interferencePhase) * factor1 * factor2;

    let probability = term1_pdg + term2_pdg + term3_pdg;

    // Ensure probability is within physical bounds [0, 1]
    probability = Math.max(0, Math.min(1, probability));

    return probability;
}

// Function to update the plot
function updatePlot() {
    // 1. Update physics parameters based on current input values
    updatePhysicsParameters();

    // 2. Read other parameters
    const deltaCp = parseFloat(dcpInput.value);
    const distance = parseFloat(distanceInput.value);
    const density = parseFloat(densityInput.value);

    // 3. Update display spans (listeners also do this, but good for consistency)
    // (Physics param spans updated in their listeners)
    dcpValueSpan.textContent = deltaCp;
    distanceValueSpan.textContent = distance;
    densityValueSpan.textContent = parseFloat(density).toFixed(1);

    // 4. Recalculate plot data
    const minEnergy = 0.2;
    const maxEnergy = 10.0;
    const numPoints = 800; // Increased number of points
    const energyStep = (maxEnergy - minEnergy) / (numPoints - 1);

    const energyValues = [];
    const probabilityValuesNu = [];
    const probabilityValuesAntiNu = [];

    for (let i = 0; i < numPoints; i++) {
        const energy = minEnergy + i * energyStep;
        energyValues.push(energy);  // Push numeric energy values for proper linear scaling

        probabilityValuesNu.push(calculateOscillationProbabilityMatter(energy, distance, density, deltaCp, false));
        probabilityValuesAntiNu.push(calculateOscillationProbabilityMatter(energy, distance, density, deltaCp, true));
    }

    // Update chart data
    oscillationChart.data.labels = energyValues;
    oscillationChart.data.datasets[0].data = probabilityValuesNu;
    oscillationChart.data.datasets[1].data = probabilityValuesAntiNu;

    // Update chart title with current parameters (using Unicode subscripts for CP)
    // Unicode subscript C: \u208C, subscript P: \u209A
    oscillationChart.options.plugins.title.text = `Oscillation Probability vs. Energy (L=${distance} km, ρ=${density.toFixed(1)} g/cm³, δCP=${deltaCp}°)`;

    oscillationChart.update('none'); // Redraw chart without animation
}

// --- Initialization ---
function initializeChart() {
    // 1. Initial update of physics parameters from default HTML values
    updatePhysicsParameters();

    // 2. Set initial span values for ALL controls
    dcpValueSpan.textContent = dcpInput.value;
    distanceValueSpan.textContent = distanceInput.value;
    densityValueSpan.textContent = parseFloat(densityInput.value).toFixed(1);
    theta12ValueSpan.textContent = parseFloat(theta12Input.value).toFixed(2);
    theta13ValueSpan.textContent = parseFloat(theta13Input.value).toFixed(2);
    theta23ValueSpan.textContent = parseFloat(theta23Input.value).toFixed(1);
    dm2_21ValueSpan.textContent = parseFloat(dm2_21Input.value).toFixed(2);
    dm2_31ValueSpan.textContent = parseFloat(dm2_31Input.value).toFixed(3);


    // 3. Calculate initial y-axis max based on default parameters
    const defaultDistance = parseFloat(distanceInput.value);
    const defaultDeltaCp = parseFloat(dcpInput.value);
    const defaultDensity = parseFloat(densityInput.value);
    let maxProb = 0.01;
    const minEnergy = 0.2;
    const maxEnergy = 10.0;
    const numPoints = 200;
    const energyStep = (maxEnergy - minEnergy) / (numPoints - 1);

    for (let i = 0; i < numPoints; i++) {
        const energy = minEnergy + i * energyStep;
        const probNu = calculateOscillationProbabilityMatter(energy, defaultDistance, defaultDensity, defaultDeltaCp, false);
        const probAntiNu = calculateOscillationProbabilityMatter(energy, defaultDistance, defaultDensity, defaultDeltaCp, true);
        maxProb = Math.max(maxProb, probNu, probAntiNu);
    }
    const yAxisMax = Math.min(1.0, Math.ceil(maxProb * 10) / 10 + 0.02); // Adjust ceiling based on typical probability ranges

    oscillationChart = new Chart(ctx, {
        type: 'line',
        data: {
            labels: [],
            datasets: [
                {
                    label: 'P(νμ → νe)', data: [],
                    borderColor: 'rgb(255, 99, 132)', backgroundColor: 'rgba(255, 99, 132, 0.1)', // Swapped to Red
                    borderWidth: 2, pointRadius: 0, tension: 0.1
                },
                {
                    label: 'P(ν̅μ → ν̅e)', data: [],
                    borderColor: 'rgb(54, 162, 235)', backgroundColor: 'rgba(54, 162, 235, 0.1)', // Swapped to Blue
                    borderWidth: 2, pointRadius: 0, tension: 0.1
                }
            ]
        },
        options: {
            responsive: true, maintainAspectRatio: true,
            scales: {
                x: {
                    type: 'linear',
                    title: {
                        display: true,
                        text: 'Neutrino Energy (E) [GeV]',
                        font: {
                            size: 16 // Increased X-axis title font size
                        }
                    },
                    min: 0,
                    max: 6
                },
                y: {
                    title: {
                        display: true,
                        text: 'Oscillation Probability',
                        font: {
                            size: 16 // Increased Y-axis title font size
                        }
                    },
                    min: 0,
                    max: 0.3
                }
            },
            plugins: {
                title: {
                    display: true,
                    text: 'Neutrino Oscillation Probability vs. Energy',
                    font: {
                        size: 18 // Increased main title font size
                    }
                },
                legend: { display: true, position: 'top' },
                tooltip: {
                    mode: 'index', intersect: false,
                    callbacks: {
                        label: function(context) {
                            let label = context.dataset.label || '';
                            if (label) label += ': ';
                            if (context.parsed.y !== null) label += `P ≈ ${context.parsed.y.toFixed(4)}`;
                            return label;
                        },
                        title: function(context) {
                            return `E = ${context[0].parsed.x} GeV`;
                        }
                    }
                },
                zoom: { // Zoom plugin configuration
                    pan: {
                        enabled: false, // Enable panning
                        mode: 'xy', // Allow panning on both axes
                        threshold: 10, // Minimum distance to initiate panning
                    },
                    zoom: {
                        wheel: {
                            enabled: false, // Disable zooming via mouse wheel
                        },
                        drag: {
                            enabled: true, // Enable drag-to-zoom (box zoom)
                            backgroundColor: 'rgba(160,160,160,0.3)', // Optional: style the drag box
                            borderColor: 'rgba(100,100,100,0.8)'
                        },
                        pinch: {
                            enabled: true // Keep pinch gesture zooming enabled
                        },
                        mode: 'xy', // Allow zooming on both axes
                    }
                }
            },
            hover: { mode: 'nearest', intersect: false },
            animation: { duration: 0 },
            datasets: { line: {} }
        }
    });

    // Add event listeners for interactivity
    dcpSlider.addEventListener('input', updatePlot);
    distanceInput.addEventListener('input', updatePlot);
    densityInput.addEventListener('input', updatePlot);

    // Initial plot
    updatePlot();

    // Ensure animations are turned off
    setTimeout(() => {
        if (oscillationChart) {
            oscillationChart.options.animation = false;
        }
    }, 100);

    // --- Synchronization Logic ---

    // Helper function for generic synchronization
    function setupSyncListeners(slider, input, valueSpan, formatFn = (v) => v) {
        slider.addEventListener('input', () => {
            const value = slider.value;
            input.value = value;
            valueSpan.textContent = formatFn(value);
            updatePlot();
        });
        input.addEventListener('input', () => {
            const value = input.value;
            const min = parseFloat(slider.min);
            const max = parseFloat(slider.max);
            const step = parseFloat(slider.step);
            let numValue = parseFloat(value);

            // Basic validation / clamping
            if (isNaN(numValue)) {
                numValue = min; // Or some default
            }
            const clampedValue = Math.max(min, Math.min(max, numValue));

            // Adjust slider to nearest step
            const steppedValue = Math.round(clampedValue / step) * step;
            // Use appropriate precision based on step for slider value setting
            const precision = (String(step).split('.')[1] || '').length;
            slider.value = steppedValue.toFixed(precision);

            // Update input field if clamped or invalid, maintaining attempted precision if possible
             const inputPrecision = (String(value).split('.')[1] || '').length;
             const displayPrecision = Math.max(precision, inputPrecision); // Use input precision for display if greater
             const finalClamped = parseFloat(clampedValue.toFixed(displayPrecision)); // Clamp with potentially higher precision

             if (input.value !== String(finalClamped) && !isNaN(finalClamped)) {
                 input.value = finalClamped;
             } else if (isNaN(finalClamped)) {
                 input.value = min; // Fallback if input was totally invalid
             }

            valueSpan.textContent = formatFn(input.value); // Update span based on potentially corrected input value
            updatePlot();
        });
    }

    // Setup listeners for all controls
    setupSyncListeners(dcpSlider, dcpInput, dcpValueSpan);
    setupSyncListeners(distanceSlider, distanceInput, distanceValueSpan);
    setupSyncListeners(densitySlider, densityInput, densityValueSpan, v => parseFloat(v).toFixed(1));
    setupSyncListeners(theta12Slider, theta12Input, theta12ValueSpan, v => parseFloat(v).toFixed(2));
    setupSyncListeners(theta13Slider, theta13Input, theta13ValueSpan, v => parseFloat(v).toFixed(2));
    setupSyncListeners(theta23Slider, theta23Input, theta23ValueSpan, v => parseFloat(v).toFixed(1));
    setupSyncListeners(dm2_21Slider, dm2_21Input, dm2_21ValueSpan, v => parseFloat(v).toFixed(2));
    setupSyncListeners(dm2_31Slider, dm2_31Input, dm2_31ValueSpan, v => parseFloat(v).toFixed(3));

    // Add listener for the mass ordering toggle
    massOrderToggle.addEventListener('change', () => {
        updatePlot(); // This will call updatePhysicsParameters which reads the toggle
    });

    // Add listener for the reset zoom button (Corrected placement and single instance)
    const resetZoomButton = document.getElementById('resetZoomButton');
    if (resetZoomButton && oscillationChart) { // Ensure button and chart exist
        resetZoomButton.addEventListener('click', () => {
            oscillationChart.resetZoom();
        });
    } // Corrected: Closing brace for the if statement

    // 5. Initial plot call
    updatePlot();
} // Corrected: End of initializeChart function

// Initialize the chart when the DOM is loaded
document.addEventListener('DOMContentLoaded', initializeChart);
// Removed duplicated block and extra closing brace
