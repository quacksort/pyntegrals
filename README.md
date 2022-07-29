# Pyntegrals
Pyntegrals is a simple library for computing a definite integral over a 2D polygonal domain.<br/>
At the moment the only function of the library is `integrate_over_polygon`, which perform the calculation.<br/>
It is possible that in the future other functions will be added.

## pyntegrals.integrals

### pyntegrals.integrals.integrate_over_polygon
    pyntegrals.integrals.integrate_over_polygon(function: Union[Callable, float], polygon: list[tuple]) -> float
 
Computes the integral of a two-variables function over a polygonal domain in R<sup>2</sup>.

Returns the definite integral of ``function(x, y)`` = f(x, y) over
``polygon`` = Σ, that is:
<div style="display: flex;justify-content: center;align-items: center;">
    <div style="display: inline-block"><span style="font-size: xx-large">∫</span><sub>Σ</sub> f(x, y) dΣ</div>
</div>

#### _Parameters_

**function** : Union[Callable, float]
<div style="margin-left: 2em;">
    A scalar or a Python function with two arguments: the first one is x and the
    second one is y in R<sup>2</sup> space.
</div>

**polygon**: list[tuple]
<div style="margin-left: 2em;">
    A polygon defined by a list of its vertices. Each vertex is
    described by a tuple containing its coordinated in the form (x, y).
    The vertices in the list must be in clockwise or
    counter-clockwise order (it makes no difference)
</div>

#### _Returns_

**r** : float
<div style="margin-left: 2em;">
    The definite integral of <code>function(x, y)</code> over <code>polygon</code>.
</div>

#### _Example_

Compute the integral of x<sup>2</sup> + y over a 50 x 40 rectangle Σ with a
vertex in (0, 0) and another in (50, 40)
<div style="display: flex;justify-content: center;align-items: center;">
    <div style="display: inline-block"><span style="font-size: xx-large">∫</span><sub>Σ</sub> (x<sup>2</sup> + y) dΣ</div>
</div>

    from pyntegrals.integrals import integrate_over_polygon
    polygon = [(0, 0), (50, 0), (50, 40), (0, 40)]
    rho = lambda x, y: x**2 + y
    mass = integrate_over_polygon(function=rho, polygon=polygon)
    >>> 1706666.666666667

For more example see examples/example.py

