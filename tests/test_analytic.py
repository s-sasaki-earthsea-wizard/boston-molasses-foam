"""Tests of the analytic reference solutions."""

from __future__ import annotations

import numpy as np
import pytest

from post.analytic import HUPPERT_ETA_N, huppert, huppert_front

AREA, G, NU = 144.0, 9.81, 55.0 / 1430.0  # Step 4-5 column, central molasses


@pytest.mark.parametrize("t", [0.0, -1.0])
def test_huppert_needs_positive_time(t: float) -> None:
    with pytest.raises(ValueError, match="t > 0"):
        huppert(np.array([0.0, 10.0]), t, AREA, G, NU)


@pytest.mark.parametrize("t", [30.0, 120.0, 600.0])
def test_huppert_profile_holds_the_area(t: float) -> None:
    x_n = huppert_front(t, AREA, G, NU)
    x = np.linspace(0.0, float(x_n), 200001)
    depth = huppert(x, t, AREA, G, NU).depth
    assert np.trapezoid(depth, x) == pytest.approx(AREA, rel=1e-4)
    assert huppert(x, t, AREA, G, NU).x_front == pytest.approx(x_n)
    assert HUPPERT_ETA_N == pytest.approx(1.41124477, rel=1e-8)


@pytest.mark.parametrize("t", [30.0, 120.0, 600.0])
def test_huppert_front_at_a_depth(t: float) -> None:
    x_d = float(huppert_front(t, AREA, G, NU, depth=0.25))
    assert x_d < float(huppert_front(t, AREA, G, NU))
    assert huppert(np.array([x_d]), t, AREA, G, NU).depth[0] == pytest.approx(0.25, rel=1e-9)
