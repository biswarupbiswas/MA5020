// Live Online Visitors Counter
function fetchOnlineCount() {
    window.WAU_r_s = function(count) {
        const onlineEl = document.getElementById('onlineCount');
        if (onlineEl) {
            const numericCount = parseInt(count.toString().replace(/,/g, ''), 10) || 1;
            onlineEl.textContent = numericCount;
        }
    };
    const script = document.createElement('script');
    script.src = 'https://whos.amung.us/pingjs/?k=ma5020mu&c=s&r=' + Math.floor(Math.random() * 10000);
    script.async = true;
    document.body.appendChild(script);
    script.onload = () => script.remove();
}

// Smooth Fast Count-Up Animation
function animateCountUp(element, target, duration = 800) {
    if (!element) return;
    const startTime = performance.now();
    function easeOutCubic(t) { return 1 - Math.pow(1 - t, 3); }

    function update(now) {
        const elapsed = now - startTime;
        const progress = Math.min(elapsed / duration, 1);
        const current = Math.floor(target * easeOutCubic(progress));
        element.textContent = current.toLocaleString('en-US');
        if (progress < 1) {
            requestAnimationFrame(update);
        } else {
            element.textContent = target.toLocaleString('en-US');
        }
    }
    requestAnimationFrame(update);
}

// Total Visitor Counter (Isolated Key: biswarupbiswas_ma5020)
function fetchTotalVisitors() {
    const totalEl = document.getElementById('totalVisitorCount');
    if (!totalEl) return;

    fetch('https://countapi.mileshilliard.com/api/v1/hit/biswarupbiswas_ma5020')
        .then(res => res.json())
        .then(data => {
            if (data && typeof data.value === 'number') {
                animateCountUp(totalEl, data.value, 800);
            }
        })
        .catch(() => {});
}

function initCounters() {
    fetchOnlineCount();
    fetchTotalVisitors();
}

if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initCounters);
} else {
    initCounters();
}

setInterval(fetchOnlineCount, 30000);
