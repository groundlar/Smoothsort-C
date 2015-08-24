Smoothsort-C
============

Implementation of Dijkstra's Smoothsort algorithm in C, for use with strings.
Current functionality is limited: the current version only works on arrays of strings
(2D arrays: arrays of pointers to arrays of chars), and only when the same case
convention is followed for every character of every string
(i.e. ALLCAPS, lowercase, aLtErNaTiNgCaSe are all fine, however CamelCase and pascalCase are not).
