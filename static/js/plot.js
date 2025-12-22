let currentData = null;
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

document.getElementById('analyze-btn').addEventListener('click', analyze);

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
        if (result.error) {
            alert('Error: ' + result.error);
        } else {
            currentData = result;
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

function displayResults(data) {
    document.getElementById('results-section').classList.remove('hidden');
    document.getElementById('display-pdb-id').textContent = data.pdb_id.toUpperCase();

    updateStats(data.phi_psi);
    generateFilters(data.phi_psi);
    renderPlot();
}

function updateStats(phiPsi) {
    const total = phiPsi.length;
    const favoured = phiPsi.filter(p => p.classification === 'favoured').length;
    const allowed = phiPsi.filter(p => p.classification === 'allowed').length;
    const outliers = phiPsi.filter(p => p.classification === 'outlier').length;

    const fPerc = total > 0 ? (favoured / total * 100).toFixed(1) : 0;
    const aPerc = total > 0 ? (allowed / total * 100).toFixed(1) : 0;
    const oPerc = total > 0 ? (outliers / total * 100).toFixed(1) : 0;

    document.getElementById('favoured-count').textContent = favoured;
    document.getElementById('favoured-percent').textContent = `(${fPerc}%)`;

    document.getElementById('allowed-count').textContent = allowed;
    document.getElementById('allowed-percent').textContent = `(${aPerc}%)`;

    document.getElementById('outlier-count').textContent = outliers;
    document.getElementById('outlier-percent').textContent = `(${oPerc}%)`;
}

function generateFilters(phiPsi) {
    const chains = [...new Set(phiPsi.map(p => p.chain))].sort();
    const types = [...new Set(phiPsi.map(p => p.rama_type))].sort();

    const chainContainer = document.getElementById('chain-filters');
    chainContainer.innerHTML = '<div class="chip active" data-filter="all">All Chains</div>';
    chains.forEach(c => {
        chainContainer.innerHTML += `<div class="chip" data-filter="${c}">${c}</div>`;
    });

    const typeContainer = document.getElementById('type-filters');
    typeContainer.innerHTML = '';
    types.forEach(t => {
        const active = currentFilters.type === t ? 'active' : '';
        typeContainer.innerHTML += `<div class="chip ${active}" data-filter="${t}">${t}</div>`;
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

function renderPlot() {
    if (!currentData) return;

    const filteredPoints = currentData.phi_psi.filter(p => {
        const chainMatch = currentFilters.chain === 'all' || p.chain === currentFilters.chain;
        const typeMatch = currentFilters.type === 'all' || p.rama_type === currentFilters.type;
        return chainMatch && typeMatch;
    });

    const traces = [];

    // 1. Background Contours
    // Use the currently selected type for contours. 
    const activeType = currentFilters.type;

    if (currentData.reference[activeType]) {
        const ref = currentData.reference[activeType];
        traces.push({
            z: ref.z,
            x: ref.phi,
            y: ref.psi,
            type: 'contour',
            showscale: false,
            contours: {
                coloring: 'none',
                start: ref.levels[0],
                end: ref.levels[1],
                size: ref.levels[1] - ref.levels[0],
                labelfont: { size: 0 }
            },
            line: {
                color: 'rgba(99, 102, 241, 0.8)', // Strengthened from 0.4
                width: 1.5 // Strengthened from 1
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
                if (p.classification === 'favoured') return '#10b981';
                if (p.classification === 'allowed') return '#f59e0b';
                return '#ef4444';
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
        font: { color: '#f8fafc', family: 'Inter' },
        xaxis: {
            title: 'Phi (Φ)',
            range: [-180, 180],
            dtick: 60,
            constrain: 'domain',
            gridcolor: 'rgba(255,255,255,0.05)',
            zerolinecolor: 'rgba(255,255,255,0.1)'
        },
        yaxis: {
            title: 'Psi (Ψ)',
            range: [-180, 180],
            dtick: 60,
            constrain: 'domain',
            gridcolor: 'rgba(255,255,255,0.05)',
            zerolinecolor: 'rgba(255,255,255,0.1)',
            scaleanchor: 'x',
            scaleratio: 1
        },
        margin: { t: 40, b: 60, l: 60, r: 40 },
        showlegend: false
    };

    const config = {
        responsive: true,
        displaylogo: false
    };

    Plotly.newPlot('rama-plot-container', traces, layout, config);
}
