(function () {
  "use strict";

  var widgets = document.querySelectorAll("[data-course-plot='joint-marginal-normal']");
  if (!widgets.length) return;

  function sequence(from, to, length) {
    var values = [];
    var step = (to - from) / (length - 1);
    for (var i = 0; i < length; i += 1) values.push(from + i * step);
    return values;
  }

  function standardNormalDensity(value) {
    return Math.exp(-0.5 * value * value) / Math.sqrt(2 * Math.PI);
  }

  function draw(widget) {
    if (widget.dataset.rendered === "true") return;

    var language = widget.dataset.lang === "en" ? "en" : "zh";
    var labels = language === "en" ? {
      title: "Joint density and its two marginal densities",
      density: "Density",
      joint: "Joint density f(x, y)",
      xMarginal: "Marginal density of X",
      yMarginal: "Marginal density of Y"
    } : {
      title: "联合密度及两个边缘密度",
      density: "密度",
      joint: "联合密度 f(x, y)",
      xMarginal: "X 的边缘密度",
      yMarginal: "Y 的边缘密度"
    };

    var x = sequence(-3.2, 3.2, 65);
    var y = sequence(-3.2, 3.2, 65);
    var xDensity = x.map(standardNormalDensity);
    var yDensity = y.map(standardNormalDensity);
    var jointDensity = y.map(function (yValue) {
      return x.map(function (xValue) {
        return standardNormalDensity(xValue) * standardNormalDensity(yValue);
      });
    });

    var traces = [
      {
        type: "surface",
        name: labels.joint,
        x: x,
        y: y,
        z: jointDensity,
        colorscale: "Viridis",
        showscale: false,
        opacity: 0.9,
        contours: {
          z: {
            show: true,
            start: 0.02,
            end: 0.16,
            size: 0.025,
            color: "rgba(35, 55, 70, 0.58)",
            project: { z: true }
          }
        },
        hovertemplate:
          "X: %{x:.2f}<br>Y: %{y:.2f}<br>" + labels.density + ": %{z:.4f}<extra>" + labels.joint + "</extra>"
      },
      {
        type: "scatter3d",
        mode: "lines",
        name: labels.xMarginal,
        x: x,
        y: x.map(function () { return -3.35; }),
        z: xDensity,
        line: { color: "#5B8FC2", width: 7 },
        hovertemplate:
          "X: %{x:.2f}<br>" + labels.density + ": %{z:.4f}<extra>" + labels.xMarginal + "</extra>"
      },
      {
        type: "scatter3d",
        mode: "lines",
        name: labels.yMarginal,
        x: y.map(function () { return 3.35; }),
        y: y,
        z: yDensity,
        line: { color: "#A43B9E", width: 7 },
        hovertemplate:
          "Y: %{y:.2f}<br>" + labels.density + ": %{z:.4f}<extra>" + labels.yMarginal + "</extra>"
      }
    ];

    var layout = {
      title: { text: labels.title, x: 0.5, xanchor: "center", font: { size: 18 } },
      autosize: true,
      margin: { l: 0, r: 0, b: 0, t: 52 },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      legend: {
        orientation: "h",
        x: 0.5,
        xanchor: "center",
        y: -0.04,
        yanchor: "top",
        bgcolor: "rgba(255,255,255,0.78)"
      },
      scene: {
        xaxis: { title: "X", range: [-3.5, 3.5], gridcolor: "#D7DEE3", zerolinecolor: "#9AA8B1" },
        yaxis: { title: "Y", range: [-3.5, 3.5], gridcolor: "#D7DEE3", zerolinecolor: "#9AA8B1" },
        zaxis: { title: labels.density, range: [0, 0.43], gridcolor: "#D7DEE3", zerolinecolor: "#9AA8B1" },
        aspectmode: "manual",
        aspectratio: { x: 1.15, y: 1.15, z: 0.9 },
        camera: { eye: { x: 1.55, y: 1.55, z: 1.05 } }
      }
    };

    var config = {
      responsive: true,
      displaylogo: false,
      scrollZoom: true,
      modeBarButtonsToRemove: ["sendDataToCloud"]
    };

    Plotly.newPlot(widget, traces, layout, config).then(function () {
      widget.dataset.rendered = "true";
      var fallback = widget.parentNode.querySelector(".plotly-fallback");
      if (fallback) fallback.hidden = true;
    });
  }

  if (typeof Plotly === "undefined") return;
  widgets.forEach(draw);
}());
