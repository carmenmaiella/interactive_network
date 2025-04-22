let fullData = {};

function createNetwork(network) {
  const width = 1000;
  const height = 700;

  const nodes = network.nodes.map(d => ({
    id: d.id,
    label: d.label,
    color: d.color
  }));

  const links = network.links.map(d => ({
    source: d.source,
    target: d.target,
    color: d.color,
    weight: d.weight,
    width: d.width,
    label: d.label || ''
  }));

  const simulation = d3.forceSimulation(nodes)
    .force("link", d3.forceLink(links).id(d => d.id).strength(0.05))
    .force("charge", d3.forceManyBody().strength(-50))
    .force("center", d3.forceCenter(width / 2, height / 2))
    .on("tick", ticked);

  simulation.alpha(1).restart();
  setTimeout(() => simulation.stop(), 3000);

  const svg = d3.create("svg")
    .attr("width", width)
    .attr("height", height)
    .attr("viewBox", [0, 0, width, height])
    .attr("style", "max-width: 100%; height: auto;");

  const containerGroup = svg.append("g");

  svg.call(d3.zoom()
    .scaleExtent([0.1, 5])
    .on("zoom", (event) => {
      containerGroup.attr("transform", event.transform);
    }));

  const link = containerGroup.selectAll("line")
    .data(links)
    .enter().append("line")
    .style("stroke", d => d.color)
    .style("stroke-width", d => d.width);

  const linkLabels = containerGroup.selectAll(".link-label")
    .data(links)
    .enter().append("text")
    .attr("class", "link-label")
    .attr("text-anchor", "middle")
    .attr("font-size", "14px")
    .attr("fill", "#555")
    .attr("font-weight", "bold")
    .style("pointer-events", "none")
    .text(d => d.label);

  const node = containerGroup.selectAll("circle")
    .data(nodes)
    .enter().append("circle")
    .attr("r", 8)
    .style("fill", d => d.color);

  node.append("title").text(d => d.label);

  node.call(d3.drag()
    .on("start", dragstarted)
    .on("drag", dragged)
    .on("end", dragended));

  function ticked() {
    link
      .attr("x1", d => d.source.x)
      .attr("y1", d => d.source.y)
      .attr("x2", d => d.target.x)
      .attr("y2", d => d.target.y);

    node
      .attr("cx", d => d.x)
      .attr("cy", d => d.y);

    linkLabels
      .attr("x", d => (d.source.x + d.target.x) / 2)
      .attr("y", d => (d.source.y + d.target.y) / 2);
  }

  function dragstarted(event) {
    simulation.alphaTarget(0.3).restart();
    event.subject.fx = event.subject.x;
    event.subject.fy = event.subject.y;
  }

  function dragged(event) {
    event.subject.fx = event.x;
    event.subject.fy = event.y;
  }

  function dragended(event) {
    event.subject.fx = null;
    event.subject.fy = null;
    simulation.alphaTarget(0);
    simulation.stop();
  }

  document.getElementById("container").innerHTML = '';
  document.getElementById("container").appendChild(svg.node());
}

document.getElementById('file-selector').addEventListener('change', (event) => {
  const file = event.target.files[0];
  if (file && file.type === 'application/json') {
    const reader = new FileReader();
    reader.onload = function(e) {
      try {
        fullData = JSON.parse(e.target.result);

        const cutoffSelect = document.getElementById('cutoff-select');
        cutoffSelect.innerHTML = '';
        Object.keys(fullData).forEach(key => {
          const option = document.createElement("option");
          option.value = key;
          option.textContent = key;
          cutoffSelect.appendChild(option);
        });
      } catch (error) {
        alert('Error parsing the JSON file. Please check the format.');
      }
    };
    reader.readAsText(file);
  } else {
    alert('Please select a valid JSON file.');
  }
});

document.getElementById("render-btn").addEventListener("click", () => {
  const selectedCutoff = document.getElementById("cutoff-select").value;
  const selectedType = document.getElementById("network-type-select").value;

  if (fullData[selectedCutoff] && fullData[selectedCutoff][selectedType]) {
    createNetwork(fullData[selectedCutoff][selectedType]);
  } else {
    alert("Network data not found for this combination.");
  }
});
