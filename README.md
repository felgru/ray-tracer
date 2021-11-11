<!--
SPDX-FileCopyrightText: 2019–2021 Felix Gruber

SPDX-License-Identifier: GPL-3.0-or-later
-->

# Ray-tracer

This is a very simple [ray tracer](https://en.wikipedia.org/wiki/Ray_tracing_(graphics)) that can render PNG files from a YAML-based scene description.

## Features

* Different geometric primitives: spheres, cubes, planes
* Procedurally generated textures: plain color, stripes, checker board
* Materials with different optical properties based on Phong shaders
* Configurable light sources
* YAML-based file format to describe scenes.

## Usage

You can build the ray tracer with Rust's `cargo` build tool, by calling
```
cargo build --release
```
(The `--release` flag is recommended as it drastically increases the speed of the renderer.)

You can play with the example scene given in `test.yaml`. Render it into
a PNG image with
```
cargo run --release test.yaml
```

## License

This program is licensed under the GPL version 3 or (at your option)
any later version.

The text of the GPL version 3 can be found in the LICENSES directory.
