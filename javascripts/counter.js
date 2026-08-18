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

// Total Visitor Counter (Actual/Original Count)
window.BusuanziCallback = function(data) {
    const totalEl = document.getElementById('totalVisitorCount');
    if (!totalEl) return;
    const rawCount = (data && data.page_pv) ? Number(data.page_pv) : ((data && data.site_pv) ? Number(data.site_pv) : 0);
    animateCountUp(totalEl, rawCount, 800);
};

function initCounters() {
    fetchOnlineCount();

    // Load live Busuanzi count with cache-busting timestamp
    const script = document.createElement('script');
    script.src = '//busuanzi.ibruce.info/busuanzi?jsonpCallback=BusuanziCallback&_=' + Date.now();
    script.async = true;
    document.body.appendChild(script);
}

if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initCounters);
} else {
    initCounters();
}

setInterval(fetchOnlineCount, 30000);
