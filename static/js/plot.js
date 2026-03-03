let currentData = null;
let cachedReferenceData = null;
let referenceRequestPromise = null;
let currentFilters = {
    chain: 'all',
    type: 'General'
};

document.getElementById('pdb-file').addEventListener('change', (e) => {
    const fileName = e.target.files[0]?.name || '';

    if (fileName) {
        // Clear PDB ID field and trigger analysis automatically
        document.getElementById('pdb-id').value = '';
        analyze();
    }
});

document.getElementById('pdb-id').addEventListener('keydown', (e) => {
    if (e.key === 'Enter') {
        analyze();
    }
});

document.getElementById('analyze-btn').addEventListener('click', analyze);
document.getElementById('write-csv-btn').addEventListener('click', downloadCSV);
document.getElementById('download-pdf-btn').addEventListener('click', downloadPDF);

async function analyze() {
    const pdbId = document.getElementById('pdb-id').value.trim();
    const pdbFile = document.getElementById('pdb-file').files[0];

    if (!pdbId && !pdbFile) {
        alert('Please provide a PDB ID or upload a file.');
        return;
    }

    const formData = new FormData();
    if (pdbId) formData.append('pdb_id', pdbId);
    if (pdbFile) formData.append('pdb_file', pdbFile);

    showLoading(true);
    try {
        const response = await fetch('/process', {
            method: 'POST',
            body: formData
        });

        const result = await response.json();
        if (!response.ok || result.error) {
            alert('Error: ' + result.error);
        } else {
            const referenceData = await ensureReferenceData();
            currentData = { ...result, reference: referenceData };
            document.getElementById('pdb-id').value = ''; // Clear on success
            displayResults(result);
        }
    } catch (err) {
        console.error(err);
        alert('An error occurred during analysis.');
    } finally {
        showLoading(false);
    }
}

function showLoading(show) {
    document.getElementById('loading').classList.toggle('hidden', !show);
}

async function ensureReferenceData() {
    if (cachedReferenceData) return cachedReferenceData;
    if (referenceRequestPromise) return referenceRequestPromise;

    referenceRequestPromise = fetch('/reference', { cache: 'force-cache' })
        .then(async (response) => {
            if (response.status === 304 && cachedReferenceData) {
                return cachedReferenceData;
            }
            if (!response.ok) {
                throw new Error('Could not load reference contour data.');
            }

            const data = await response.json();
            cachedReferenceData = data;
            return data;
        })
        .finally(() => {
            referenceRequestPromise = null;
        });

    return referenceRequestPromise;
}

function displayResults(data) {
    // Reset filters to defaults for new structure
    currentFilters.chain = 'all';
    currentFilters.type = getPreferredType(data.phi_psi);

    document.getElementById('results-section').classList.remove('hidden');
    generateFilters(data.phi_psi);
    renderPlot();
}

function getPreferredType(phiPsi) {
    const types = [...new Set(phiPsi.map(p => p.rama_type))].sort();
    if (types.includes('General')) return 'General';
    return types[0] || 'General';
}

function updateStats(phiPsi, prefix = '') {
    const validPoints = phiPsi.filter(p => p.phi !== null && p.psi !== null);
    const total = validPoints.length;
    const favoured = validPoints.filter(p => p.classification === 'favoured').length;
    const allowed = validPoints.filter(p => p.classification === 'allowed').length;
    const outliers = validPoints.filter(p => p.classification === 'outlier').length;

    const fPerc = total > 0 ? (favoured / total * 100).toFixed(1) : 0;
    const aPerc = total > 0 ? (allowed / total * 100).toFixed(1) : 0;
    const oPerc = total > 0 ? (outliers / total * 100).toFixed(1) : 0;

    document.getElementById(`${prefix}favoured-count`).textContent = favoured;
    document.getElementById(`${prefix}favoured-percent`).textContent = `(${fPerc}%)`;

    document.getElementById(`${prefix}allowed-count`).textContent = allowed;
    document.getElementById(`${prefix}allowed-percent`).textContent = `(${aPerc}%)`;

    document.getElementById(`${prefix}outlier-count`).textContent = outliers;
    document.getElementById(`${prefix}outlier-percent`).textContent = `(${oPerc}%)`;
}

function generateFilters(phiPsi) {
    const chains = [...new Set(phiPsi.map(p => p.chain))].sort();
    const types = [...new Set(phiPsi.map(p => p.rama_type))].sort();
    if (!types.includes(currentFilters.type)) {
        currentFilters.type = getPreferredType(phiPsi);
    }

    const chainContainer = document.getElementById('chain-filters');
    chainContainer.replaceChildren();
    const allChainsChip = document.createElement('div');
    allChainsChip.className = 'chip active';
    allChainsChip.dataset.filter = 'all';
    allChainsChip.textContent = 'All Chains';
    chainContainer.appendChild(allChainsChip);
    chains.forEach(c => {
        const chip = document.createElement('div');
        chip.className = 'chip';
        chip.dataset.filter = c;
        chip.textContent = c;
        chainContainer.appendChild(chip);
    });

    const typeContainer = document.getElementById('type-filters');
    typeContainer.replaceChildren();
    types.forEach(t => {
        const active = currentFilters.type === t ? 'active' : '';
        const chip = document.createElement('div');
        chip.className = `chip ${active}`.trim();
        chip.dataset.filter = t;
        chip.textContent = t;
        typeContainer.appendChild(chip);
    });

    // Add event listeners
    document.querySelectorAll('#chain-filters .chip').forEach(chip => {
        chip.addEventListener('click', () => {
            document.querySelectorAll('#chain-filters .chip').forEach(c => c.classList.remove('active'));
            chip.classList.add('active');
            currentFilters.chain = chip.dataset.filter;
            renderPlot();
        });
    });

    document.querySelectorAll('#type-filters .chip').forEach(chip => {
        chip.addEventListener('click', () => {
            document.querySelectorAll('#type-filters .chip').forEach(c => c.classList.remove('active'));
            chip.classList.add('active');
            currentFilters.type = chip.dataset.filter;
            renderPlot();
        });
    });
}

function updateSummaryHeader() {
    if (!currentData) return;

    document.getElementById('display-pdb-id').textContent = currentData.pdb_id.toUpperCase();

    const parts = [];
    if (currentFilters.chain !== 'all') {
        parts.push(`Chain ${currentFilters.chain}`);
    }
    if (currentFilters.type !== 'all') {
        parts.push(currentFilters.type);
    }

    const filterInfo = parts.length > 0 ? parts.join(' • ') : 'All Residues';
    document.getElementById('display-filter-info').textContent = filterInfo;
}

function renderPlot() {
    if (!currentData) return;

    const filteredPoints = currentData.phi_psi.filter(p => {
        if (p.phi === null || p.psi === null) return false;
        const chainMatch = currentFilters.chain === 'all' || p.chain === currentFilters.chain;
        const typeMatch = currentFilters.type === 'all' || p.rama_type === currentFilters.type;
        return chainMatch && typeMatch;
    });

    updateStats(filteredPoints, '');
    updateStats(currentData.phi_psi, 'overall-');
    updateSummaryHeader();

    const traces = [];

    // 1. Background Contours
    // Use the currently selected type for contours. 
    const activeType = currentFilters.type;

    if (currentData.reference[activeType]) {
        const ref = currentData.reference[activeType];

        // 1a. Outer Contour (Allowed Region) - Purple/Magenta
        traces.push({
            z: ref.z,
            x: ref.phi,
            y: ref.psi,
            type: 'contour',
            showscale: false,
            contours: {
                coloring: 'none',
                start: ref.levels[0],
                end: ref.levels[0],
                size: 1,
                labelfont: { size: 0 }
            },
            line: {
                color: '#7e48dbff', // Purple
                width: 1.5
            },
            hoverinfo: 'skip'
        });

        // 1b. Inner Contour (Favoured Region) - Blue
        traces.push({
            z: ref.z,
            x: ref.phi,
            y: ref.psi,
            type: 'contour',
            showscale: false,
            contours: {
                coloring: 'none',
                start: ref.levels[1],
                end: ref.levels[1],
                size: 1,
                labelfont: { size: 0 }
            },
            line: {
                color: '#3036e7ff', // Blue
                width: 1.5
            },
            hoverinfo: 'skip'
        });
    }

    // 2. Scatter Points
    const scatterTrace = {
        x: filteredPoints.map(p => p.phi),
        y: filteredPoints.map(p => p.psi),
        mode: 'markers',
        type: 'scatter',
        marker: {
            size: 8,
            color: filteredPoints.map(p => {
                if (p.classification === 'favoured') return '#34d399'; // Softer Emerald
                if (p.classification === 'allowed') return '#fbbf24'; // Softer Amber
                return '#f87171'; // Softer Red
            }),
            symbol: 'circle-open',
            line: { width: 1.5 }
        },
        text: filteredPoints.map(p => `${p.resName} ${p.resSeq} (Chain ${p.chain})<br>Phi: ${p.phi.toFixed(2)}, Psi: ${p.psi.toFixed(2)}<br>Type: ${p.rama_type}<br>Category: ${p.classification}`),
        hoverinfo: 'text'
    };
    traces.push(scatterTrace);

    const layout = {
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: { color: '#1e293b', family: 'Inter' },
        xaxis: {
            title: 'Phi (Φ)',
            range: [-180, 180],
            dtick: 60,
            constrain: 'domain',
            gridcolor: 'rgba(0,0,0,0.05)',
            zerolinecolor: 'rgba(0,0,0,0.1)'
        },
        yaxis: {
            title: 'Psi (Ψ)',
            range: [-180, 180],
            dtick: 60,
            constrain: 'domain',
            gridcolor: 'rgba(0,0,0,0.05)',
            zerolinecolor: 'rgba(0,0,0,0.1)',
            scaleanchor: 'x',
            scaleratio: 1
        },
        margin: { t: 40, b: 60, l: 60, r: 40 },
        showlegend: false
    };

    const config = {
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d']
    };

    Plotly.newPlot('rama-plot-container', traces, layout, config);
}

function downloadCSV() {
    if (!currentData || !currentData.result_id) {
        alert("No analysis result available yet.");
        return;
    }
    window.location.href = `/download/csv/${currentData.result_id}`;
}

function downloadPDF() {
    if (!currentData || !currentData.result_id) {
        alert("No analysis result available yet.");
        return;
    }
    window.location.href = `/download/pdf/${currentData.result_id}`;
}
