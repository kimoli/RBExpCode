function plotMedianBoxplot(vector, center, halfwidth, color)

tempquants = quantile(vector,3);
plot([center-halfwidth center+halfwidth], [tempquants(2), tempquants(2)], 'Color', color)
plot([center-halfwidth center+halfwidth], [tempquants(1), tempquants(1)], 'Color', color)
plot([center-halfwidth center+halfwidth], [tempquants(3), tempquants(3)], 'Color', color)
plot([center-halfwidth center-halfwidth], [tempquants(1), tempquants(3)], 'Color', color)
plot([center+halfwidth center+halfwidth], [tempquants(1), tempquants(3)], 'Color', color)
plot([center center], [tempquants(1), min(vector)], 'Color', color)
plot([center center], [max(vector), tempquants(3)], 'Color', color)

end