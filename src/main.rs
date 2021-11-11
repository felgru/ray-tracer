// SPDX-FileCopyrightText: 2019–2021 Felix Gruber
//
// SPDX-License-Identifier: GPL-3.0-or-later

use std::env;
use std::fs;
use std::path::Path;

use ray_tracer::yaml::load_world_and_cameras_from_str;

fn main() {
    if let Some(path) = env::args().nth(1) {
        let path = Path::new(&path);
        render_scene(&path);
    } else {
        println!("Need to provide filename of scene description as argument");
    }
}

fn render_scene(path: &Path) {
    let file = fs::read_to_string(path).unwrap();
    let (world, cameras) = load_world_and_cameras_from_str(&file);

    let canvas = cameras[0].render(&world);
    let render_file = Path::new(path.file_name().unwrap())
                      .with_extension("png");
    canvas.write_image(&render_file).unwrap();
}
