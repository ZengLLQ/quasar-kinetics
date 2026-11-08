package engine

import (
	"math"

	"github.com/ZengLLQ/quasar-kinetics/pkg/particles"
)

// Tile represents a contiguous spatial partition of the simulation domain.
type Tile struct {
	ID       int
	XMin     float64
	XMax     float64
	Indices  []int // Indices into the species particle slice
}

// TilePartition manages domain decomposition into spatial tiles for cache-friendly access.
type TilePartition struct {
	NumTiles  int
	TileWidth float64
	Lx        float64
	X0        float64
	Tiles     []Tile
}

// NewTilePartition creates a 1D spatial tile partitioning.
func NewTilePartition(numTiles int, x0, lx float64) *TilePartition {
	if numTiles <= 0 {
		numTiles = 4
	}
	tileWidth := lx / float64(numTiles)
	tiles := make([]Tile, numTiles)
	for i := 0; i < numTiles; i++ {
		tiles[i] = Tile{
			ID:      i,
			XMin:    x0 + float64(i)*tileWidth,
			XMax:    x0 + float64(i+1)*tileWidth,
			Indices: make([]int, 0, 1024),
		}
	}
	return &TilePartition{
		NumTiles:  numTiles,
		TileWidth: tileWidth,
		Lx:        lx,
		X0:        x0,
		Tiles:     tiles,
	}
}

// BinSpecies assigns each particle of a species to its corresponding spatial tile.
func (tp *TilePartition) BinSpecies(sp *particles.Species) {
	for i := range tp.Tiles {
		tp.Tiles[i].Indices = tp.Tiles[i].Indices[:0]
	}

	invWidth := 1.0 / tp.TileWidth
	numTiles := tp.NumTiles

	for pIdx := range sp.Particles {
		x := sp.Particles[pIdx].X
		// Normalize to [0, Lx)
		xRel := math.Mod(x-tp.X0, tp.Lx)
		if xRel < 0 {
			xRel += tp.Lx
		}

		tIdx := int(xRel * invWidth)
		if tIdx >= numTiles {
			tIdx = numTiles - 1
		} else if tIdx < 0 {
			tIdx = 0
		}

		tp.Tiles[tIdx].Indices = append(tp.Tiles[tIdx].Indices, pIdx)
	}
}
