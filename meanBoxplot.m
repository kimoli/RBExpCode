function meanBoxplot(data,center,width,color)
hold on
semval = std(data)./sqrt(length(data));
plot([center-(width/2) center+(width/2)], [mean(data), mean(data)], 'Color', color)
plot([center-(width/2) center+(width/2)], [mean(data)+semval, mean(data)+semval], 'Color', color)
plot([center-(width/2) center+(width/2)], [mean(data)-semval, mean(data)-semval], 'Color', color)
plot([center-(width/2) center-(width/2)], [mean(data)-semval, mean(data)+semval], 'Color', color)
plot([center+(width/2) center+(width/2)], [mean(data)-semval, mean(data)+semval], 'Color', color)
plot([center center], [mean(data)-semval, min(data)], 'Color', color)
plot([center center], [mean(data)+semval, max(data)], 'Color', color)
scatter(ones(length(data),1)*center, data, 15, 'MarkerEdgeColor', color)

end